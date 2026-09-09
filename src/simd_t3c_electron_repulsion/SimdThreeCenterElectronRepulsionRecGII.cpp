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


#include "SimdThreeCenterElectronRepulsionRecGII.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
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

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_gii_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_gii_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 151046, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1521 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 151046, 110403, 8572, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 16,
                                                             ncols, fj, mu, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2724, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2727, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2730, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2733, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2736, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2739, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2742, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2745, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2748, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2751, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2754, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2757, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2760, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2763, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2766, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2769, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2778, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2787, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2796, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2805, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2814, 3, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2823, 3, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2832, 3, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2841, 3, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2850, 3, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2859, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2868, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2877, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2886, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2895, 3, 30, 84,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2913, 3, 33, 90,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2931, 3, 36, 96,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2949, 3, 39, 102,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2967, 3, 42, 108,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2985, 3, 45, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3003, 3, 48, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3021, 3, 51, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3039, 3, 54, 132,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3057, 3, 57, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3075, 3, 60, 144,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3093, 3, 63, 150,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3111, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3129, 3, 84, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3159, 3, 90, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3189, 3, 96, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3219, 3, 102, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3249, 3, 108, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3279, 3, 114, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3309, 3, 120, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3339, 3, 126, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3369, 3, 132, 262,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3399, 3, 138, 272,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3429, 3, 144, 282,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3459, 3, 150, 292,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3489, 3, 182, 332,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3534, 3, 192, 347,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3579, 3, 202, 362,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3624, 3, 212, 377,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3669, 3, 222, 392,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3714, 3, 232, 407,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3759, 3, 242, 422,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3804, 3, 252, 437,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3849, 3, 262, 452,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3894, 3, 272, 467,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3939, 3, 282, 482,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3984, 3, 332, 539,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4047, 3, 347, 560,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4110, 3, 362, 581,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4173, 3, 377, 602,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4236, 3, 392, 623,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4299, 3, 407, 644,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4362, 3, 422, 665,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4425, 3, 437, 686,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4488, 3, 452, 707,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4551, 3, 467, 728,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4614, 3, 539, 805,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4698, 3, 560, 833,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4782, 3, 581, 861,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4866, 3, 602, 889,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4950, 3, 623, 917,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5034, 3, 644, 945,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5118, 3, 665, 973,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5202, 3, 686,
                                                                       1001, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5286, 3, 707,
                                                                       1029, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5370, 3, 805,
                                                                       1129, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5478, 3, 833,
                                                                       1165, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5586, 3, 861,
                                                                       1201, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5694, 3, 889,
                                                                       1237, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5802, 3, 917,
                                                                       1273, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5910, 3, 945,
                                                                       1309, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6018, 3, 973,
                                                                       1345, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6126, 3, 1001,
                                                                       1381, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6234, 3, 1129,
                                                                       1507, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6369, 3, 1165,
                                                                       1552, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6504, 3, 1201,
                                                                       1597, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6639, 3, 1237,
                                                                       1642, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6774, 3, 1273,
                                                                       1687, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6909, 3, 1309,
                                                                       1732, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7044, 3, 1345,
                                                                       1777, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7179, 3, 1507,
                                                                       1932, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7344, 3, 1552,
                                                                       1987, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7509, 3, 1597,
                                                                       2042, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7674, 3, 1642,
                                                                       2097, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7839, 3, 1687,
                                                                       2152, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8004, 3, 1732,
                                                                       2207, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8169, 3, 1932,
                                                                       2394, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8367, 3, 1987,
                                                                       2460, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8565, 3, 2042,
                                                                       2526, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8763, 3, 2097,
                                                                       2592, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8961, 3, 2152,
                                                                       2658, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9159, 3, 7, 8,
                                                                       2724, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9165, 3, 8, 9,
                                                                       2727, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9171, 3, 9, 10,
                                                                       2730, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9177, 3, 10, 11,
                                                                       2733, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9183, 3, 11, 12,
                                                                       2736, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9189, 3, 12, 13,
                                                                       2739, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9195, 3, 13, 14,
                                                                       2742, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9201, 3, 14, 15,
                                                                       2745, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9207, 3, 15, 16,
                                                                       2748, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9213, 3, 16, 17,
                                                                       2751, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9219, 3, 17, 18,
                                                                       2754, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9225, 3, 18, 19,
                                                                       2757, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9231, 3, 19, 20,
                                                                       2760, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9237, 3, 20, 21,
                                                                       2763, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9243, 3, 21, 22,
                                                                       2766, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9249, 0, 3, 9159,
                                                                       2724, 9165, 2769, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9267, 0, 3, 9165,
                                                                       2727, 9171, 2778, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9285, 0, 3, 9171,
                                                                       2730, 9177, 2787, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9303, 0, 3, 9177,
                                                                       2733, 9183, 2796, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9321, 0, 3, 9183,
                                                                       2736, 9189, 2805, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9339, 0, 3, 9189,
                                                                       2739, 9195, 2814, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9357, 0, 3, 9195,
                                                                       2742, 9201, 2823, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9375, 0, 3, 9201,
                                                                       2745, 9207, 2832, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9393, 0, 3, 9207,
                                                                       2748, 9213, 2841, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9411, 0, 3, 9213,
                                                                       2751, 9219, 2850, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9429, 0, 3, 9219,
                                                                       2754, 9225, 2859, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9447, 0, 3, 9225,
                                                                       2757, 9231, 2868, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9465, 0, 3, 9231,
                                                                       2760, 9237, 2877, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9483, 0, 3, 9237,
                                                                       2763, 9243, 2886, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9501, 0, 3, 9249,
                                                                       2769, 9267, 72, 78, 2895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9537, 0, 3, 9267,
                                                                       2778, 9285, 78, 84, 2913,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9573, 0, 3, 9285,
                                                                       2787, 9303, 84, 90, 2931,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9609, 0, 3, 9303,
                                                                       2796, 9321, 90, 96, 2949,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9645, 0, 3, 9321,
                                                                       2805, 9339, 96, 102, 2967,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9681, 0, 3, 9339,
                                                                       2814, 9357, 102, 108,
                                                                       2985, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9717, 0, 3, 9357,
                                                                       2823, 9375, 108, 114,
                                                                       3003, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9753, 0, 3, 9375,
                                                                       2832, 9393, 114, 120,
                                                                       3021, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9789, 0, 3, 9393,
                                                                       2841, 9411, 120, 126,
                                                                       3039, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9825, 0, 3, 9411,
                                                                       2850, 9429, 126, 132,
                                                                       3057, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9861, 0, 3, 9429,
                                                                       2859, 9447, 132, 138,
                                                                       3075, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9897, 0, 3, 9447,
                                                                       2868, 9465, 138, 144,
                                                                       3093, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9933, 0, 3, 9465,
                                                                       2877, 9483, 144, 150,
                                                                       3111, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9969, 0, 3, 9501,
                                                                       2895, 9537, 162, 172,
                                                                       3129, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10029, 0, 3, 9537,
                                                                       2913, 9573, 172, 182,
                                                                       3159, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10089, 0, 3, 9573,
                                                                       2931, 9609, 182, 192,
                                                                       3189, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10149, 0, 3, 9609,
                                                                       2949, 9645, 192, 202,
                                                                       3219, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10209, 0, 3, 9645,
                                                                       2967, 9681, 202, 212,
                                                                       3249, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10269, 0, 3, 9681,
                                                                       2985, 9717, 212, 222,
                                                                       3279, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10329, 0, 3, 9717,
                                                                       3003, 9753, 222, 232,
                                                                       3309, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10389, 0, 3, 9753,
                                                                       3021, 9789, 232, 242,
                                                                       3339, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10449, 0, 3, 9789,
                                                                       3039, 9825, 242, 252,
                                                                       3369, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10509, 0, 3, 9825,
                                                                       3057, 9861, 252, 262,
                                                                       3399, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10569, 0, 3, 9861,
                                                                       3075, 9897, 262, 272,
                                                                       3429, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10629, 0, 3, 9897,
                                                                       3093, 9933, 272, 282,
                                                                       3459, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10689, 0, 3, 9969,
                                                                       3129, 10029, 302, 317,
                                                                       3489, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10779, 0, 3,
                                                                       10029, 3159, 10089, 317,
                                                                       332, 3534, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10869, 0, 3,
                                                                       10089, 3189, 10149, 332,
                                                                       347, 3579, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10959, 0, 3,
                                                                       10149, 3219, 10209, 347,
                                                                       362, 3624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11049, 0, 3,
                                                                       10209, 3249, 10269, 362,
                                                                       377, 3669, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11139, 0, 3,
                                                                       10269, 3279, 10329, 377,
                                                                       392, 3714, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11229, 0, 3,
                                                                       10329, 3309, 10389, 392,
                                                                       407, 3759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11319, 0, 3,
                                                                       10389, 3339, 10449, 407,
                                                                       422, 3804, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11409, 0, 3,
                                                                       10449, 3369, 10509, 422,
                                                                       437, 3849, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11499, 0, 3,
                                                                       10509, 3399, 10569, 437,
                                                                       452, 3894, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11589, 0, 3,
                                                                       10569, 3429, 10629, 452,
                                                                       467, 3939, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11679, 0, 3,
                                                                       10689, 3489, 10779, 497,
                                                                       518, 3984, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11805, 0, 3,
                                                                       10779, 3534, 10869, 518,
                                                                       539, 4047, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11931, 0, 3,
                                                                       10869, 3579, 10959, 539,
                                                                       560, 4110, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12057, 0, 3,
                                                                       10959, 3624, 11049, 560,
                                                                       581, 4173, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12183, 0, 3,
                                                                       11049, 3669, 11139, 581,
                                                                       602, 4236, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12309, 0, 3,
                                                                       11139, 3714, 11229, 602,
                                                                       623, 4299, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12435, 0, 3,
                                                                       11229, 3759, 11319, 623,
                                                                       644, 4362, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12561, 0, 3,
                                                                       11319, 3804, 11409, 644,
                                                                       665, 4425, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12687, 0, 3,
                                                                       11409, 3849, 11499, 665,
                                                                       686, 4488, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12813, 0, 3,
                                                                       11499, 3894, 11589, 686,
                                                                       707, 4551, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12939, 0, 3,
                                                                       11679, 3984, 11805, 749,
                                                                       777, 4614, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13107, 0, 3,
                                                                       11805, 4047, 11931, 777,
                                                                       805, 4698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13275, 0, 3,
                                                                       11931, 4110, 12057, 805,
                                                                       833, 4782, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13443, 0, 3,
                                                                       12057, 4173, 12183, 833,
                                                                       861, 4866, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13611, 0, 3,
                                                                       12183, 4236, 12309, 861,
                                                                       889, 4950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13779, 0, 3,
                                                                       12309, 4299, 12435, 889,
                                                                       917, 5034, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13947, 0, 3,
                                                                       12435, 4362, 12561, 917,
                                                                       945, 5118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14115, 0, 3,
                                                                       12561, 4425, 12687, 945,
                                                                       973, 5202, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14283, 0, 3,
                                                                       12687, 4488, 12813, 973,
                                                                       1001, 5286, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14451, 0, 3,
                                                                       12939, 4614, 13107, 1057,
                                                                       1093, 5370, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14667, 0, 3,
                                                                       13107, 4698, 13275, 1093,
                                                                       1129, 5478, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14883, 0, 3,
                                                                       13275, 4782, 13443, 1129,
                                                                       1165, 5586, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15099, 0, 3,
                                                                       13443, 4866, 13611, 1165,
                                                                       1201, 5694, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15315, 0, 3,
                                                                       13611, 4950, 13779, 1201,
                                                                       1237, 5802, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15531, 0, 3,
                                                                       13779, 5034, 13947, 1237,
                                                                       1273, 5910, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15747, 0, 3,
                                                                       13947, 5118, 14115, 1273,
                                                                       1309, 6018, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15963, 0, 3,
                                                                       14115, 5202, 14283, 1309,
                                                                       1345, 6126, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16179, 0, 3,
                                                                       14451, 5370, 14667, 1417,
                                                                       1462, 6234, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16449, 0, 3,
                                                                       14667, 5478, 14883, 1462,
                                                                       1507, 6369, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16719, 0, 3,
                                                                       14883, 5586, 15099, 1507,
                                                                       1552, 6504, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16989, 0, 3,
                                                                       15099, 5694, 15315, 1552,
                                                                       1597, 6639, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17259, 0, 3,
                                                                       15315, 5802, 15531, 1597,
                                                                       1642, 6774, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17529, 0, 3,
                                                                       15531, 5910, 15747, 1642,
                                                                       1687, 6909, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17799, 0, 3,
                                                                       15747, 6018, 15963, 1687,
                                                                       1732, 7044, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 18069, 0, 3,
                                                                       16179, 6234, 16449, 1822,
                                                                       1877, 7179, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 18399, 0, 3,
                                                                       16449, 6369, 16719, 1877,
                                                                       1932, 7344, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 18729, 0, 3,
                                                                       16719, 6504, 16989, 1932,
                                                                       1987, 7509, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19059, 0, 3,
                                                                       16989, 6639, 17259, 1987,
                                                                       2042, 7674, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19389, 0, 3,
                                                                       17259, 6774, 17529, 2042,
                                                                       2097, 7839, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19719, 0, 3,
                                                                       17529, 6909, 17799, 2097,
                                                                       2152, 8004, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 20049, 0, 3,
                                                                       18069, 7179, 18399, 2262,
                                                                       2328, 8169, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 20445, 0, 3,
                                                                       18399, 7344, 18729, 2328,
                                                                       2394, 8367, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 20841, 0, 3,
                                                                       18729, 7509, 19059, 2394,
                                                                       2460, 8565, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 21237, 0, 3,
                                                                       19059, 7674, 19389, 2460,
                                                                       2526, 8763, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 21633, 0, 3,
                                                                       19389, 7839, 19719, 2526,
                                                                       2592, 8961, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22029, 3, 2724,
                                                                       2727, 9171, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22039, 3, 2727,
                                                                       2730, 9177, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22049, 3, 2730,
                                                                       2733, 9183, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22059, 3, 2733,
                                                                       2736, 9189, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22069, 3, 2736,
                                                                       2739, 9195, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22079, 3, 2739,
                                                                       2742, 9201, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22089, 3, 2742,
                                                                       2745, 9207, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22099, 3, 2745,
                                                                       2748, 9213, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22109, 3, 2748,
                                                                       2751, 9219, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22119, 3, 2751,
                                                                       2754, 9225, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22129, 3, 2754,
                                                                       2757, 9231, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22139, 3, 2757,
                                                                       2760, 9237, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22149, 3, 2760,
                                                                       2763, 9243, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22159, 0, 3,
                                                                       22029, 9171, 22039, 9285,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22189, 0, 3,
                                                                       22039, 9177, 22049, 9303,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22219, 0, 3,
                                                                       22049, 9183, 22059, 9321,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22249, 0, 3,
                                                                       22059, 9189, 22069, 9339,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22279, 0, 3,
                                                                       22069, 9195, 22079, 9357,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22309, 0, 3,
                                                                       22079, 9201, 22089, 9375,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22339, 0, 3,
                                                                       22089, 9207, 22099, 9393,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22369, 0, 3,
                                                                       22099, 9213, 22109, 9411,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22399, 0, 3,
                                                                       22109, 9219, 22119, 9429,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22429, 0, 3,
                                                                       22119, 9225, 22129, 9447,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22459, 0, 3,
                                                                       22129, 9231, 22139, 9465,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22489, 0, 3,
                                                                       22139, 9237, 22149, 9483,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22519, 0, 3,
                                                                       22159, 9285, 22189, 2895,
                                                                       2913, 9573, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22579, 0, 3,
                                                                       22189, 9303, 22219, 2913,
                                                                       2931, 9609, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22639, 0, 3,
                                                                       22219, 9321, 22249, 2931,
                                                                       2949, 9645, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22699, 0, 3,
                                                                       22249, 9339, 22279, 2949,
                                                                       2967, 9681, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22759, 0, 3,
                                                                       22279, 9357, 22309, 2967,
                                                                       2985, 9717, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22819, 0, 3,
                                                                       22309, 9375, 22339, 2985,
                                                                       3003, 9753, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22879, 0, 3,
                                                                       22339, 9393, 22369, 3003,
                                                                       3021, 9789, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22939, 0, 3,
                                                                       22369, 9411, 22399, 3021,
                                                                       3039, 9825, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22999, 0, 3,
                                                                       22399, 9429, 22429, 3039,
                                                                       3057, 9861, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23059, 0, 3,
                                                                       22429, 9447, 22459, 3057,
                                                                       3075, 9897, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23119, 0, 3,
                                                                       22459, 9465, 22489, 3075,
                                                                       3093, 9933, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23179, 0, 3,
                                                                       22519, 9573, 22579, 3129,
                                                                       3159, 10089, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23279, 0, 3,
                                                                       22579, 9609, 22639, 3159,
                                                                       3189, 10149, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23379, 0, 3,
                                                                       22639, 9645, 22699, 3189,
                                                                       3219, 10209, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23479, 0, 3,
                                                                       22699, 9681, 22759, 3219,
                                                                       3249, 10269, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23579, 0, 3,
                                                                       22759, 9717, 22819, 3249,
                                                                       3279, 10329, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23679, 0, 3,
                                                                       22819, 9753, 22879, 3279,
                                                                       3309, 10389, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23779, 0, 3,
                                                                       22879, 9789, 22939, 3309,
                                                                       3339, 10449, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23879, 0, 3,
                                                                       22939, 9825, 22999, 3339,
                                                                       3369, 10509, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23979, 0, 3,
                                                                       22999, 9861, 23059, 3369,
                                                                       3399, 10569, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24079, 0, 3,
                                                                       23059, 9897, 23119, 3399,
                                                                       3429, 10629, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24179, 0, 3,
                                                                       23179, 10089, 23279, 3489,
                                                                       3534, 10869, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24329, 0, 3,
                                                                       23279, 10149, 23379, 3534,
                                                                       3579, 10959, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24479, 0, 3,
                                                                       23379, 10209, 23479, 3579,
                                                                       3624, 11049, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24629, 0, 3,
                                                                       23479, 10269, 23579, 3624,
                                                                       3669, 11139, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24779, 0, 3,
                                                                       23579, 10329, 23679, 3669,
                                                                       3714, 11229, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24929, 0, 3,
                                                                       23679, 10389, 23779, 3714,
                                                                       3759, 11319, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25079, 0, 3,
                                                                       23779, 10449, 23879, 3759,
                                                                       3804, 11409, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25229, 0, 3,
                                                                       23879, 10509, 23979, 3804,
                                                                       3849, 11499, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25379, 0, 3,
                                                                       23979, 10569, 24079, 3849,
                                                                       3894, 11589, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25529, 0, 3,
                                                                       24179, 10869, 24329, 3984,
                                                                       4047, 11931, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25739, 0, 3,
                                                                       24329, 10959, 24479, 4047,
                                                                       4110, 12057, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25949, 0, 3,
                                                                       24479, 11049, 24629, 4110,
                                                                       4173, 12183, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26159, 0, 3,
                                                                       24629, 11139, 24779, 4173,
                                                                       4236, 12309, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26369, 0, 3,
                                                                       24779, 11229, 24929, 4236,
                                                                       4299, 12435, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26579, 0, 3,
                                                                       24929, 11319, 25079, 4299,
                                                                       4362, 12561, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26789, 0, 3,
                                                                       25079, 11409, 25229, 4362,
                                                                       4425, 12687, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26999, 0, 3,
                                                                       25229, 11499, 25379, 4425,
                                                                       4488, 12813, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27209, 0, 3,
                                                                       25529, 11931, 25739, 4614,
                                                                       4698, 13275, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27489, 0, 3,
                                                                       25739, 12057, 25949, 4698,
                                                                       4782, 13443, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27769, 0, 3,
                                                                       25949, 12183, 26159, 4782,
                                                                       4866, 13611, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28049, 0, 3,
                                                                       26159, 12309, 26369, 4866,
                                                                       4950, 13779, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28329, 0, 3,
                                                                       26369, 12435, 26579, 4950,
                                                                       5034, 13947, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28609, 0, 3,
                                                                       26579, 12561, 26789, 5034,
                                                                       5118, 14115, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28889, 0, 3,
                                                                       26789, 12687, 26999, 5118,
                                                                       5202, 14283, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29169, 0, 3,
                                                                       27209, 13275, 27489, 5370,
                                                                       5478, 14883, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29529, 0, 3,
                                                                       27489, 13443, 27769, 5478,
                                                                       5586, 15099, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29889, 0, 3,
                                                                       27769, 13611, 28049, 5586,
                                                                       5694, 15315, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30249, 0, 3,
                                                                       28049, 13779, 28329, 5694,
                                                                       5802, 15531, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30609, 0, 3,
                                                                       28329, 13947, 28609, 5802,
                                                                       5910, 15747, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30969, 0, 3,
                                                                       28609, 14115, 28889, 5910,
                                                                       6018, 15963, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 31329, 0, 3,
                                                                       29169, 14883, 29529, 6234,
                                                                       6369, 16719, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 31779, 0, 3,
                                                                       29529, 15099, 29889, 6369,
                                                                       6504, 16989, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 32229, 0, 3,
                                                                       29889, 15315, 30249, 6504,
                                                                       6639, 17259, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 32679, 0, 3,
                                                                       30249, 15531, 30609, 6639,
                                                                       6774, 17529, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 33129, 0, 3,
                                                                       30609, 15747, 30969, 6774,
                                                                       6909, 17799, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 33579, 0, 3,
                                                                       31329, 16719, 31779, 7179,
                                                                       7344, 18729, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 34129, 0, 3,
                                                                       31779, 16989, 32229, 7344,
                                                                       7509, 19059, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 34679, 0, 3,
                                                                       32229, 17259, 32679, 7509,
                                                                       7674, 19389, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 35229, 0, 3,
                                                                       32679, 17529, 33129, 7674,
                                                                       7839, 19719, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 35779, 0, 3,
                                                                       33579, 18729, 34129, 8169,
                                                                       8367, 20841, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 36439, 0, 3,
                                                                       34129, 19059, 34679, 8367,
                                                                       8565, 21237, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 37099, 0, 3,
                                                                       34679, 19389, 35229, 8565,
                                                                       8763, 21633, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37759, 3, 9159,
                                                                       9165, 22029, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37774, 3, 9165,
                                                                       9171, 22039, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37789, 3, 9171,
                                                                       9177, 22049, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37804, 3, 9177,
                                                                       9183, 22059, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37819, 3, 9183,
                                                                       9189, 22069, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37834, 3, 9189,
                                                                       9195, 22079, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37849, 3, 9195,
                                                                       9201, 22089, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37864, 3, 9201,
                                                                       9207, 22099, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37879, 3, 9207,
                                                                       9213, 22109, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37894, 3, 9213,
                                                                       9219, 22119, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37909, 3, 9219,
                                                                       9225, 22129, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37924, 3, 9225,
                                                                       9231, 22139, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37939, 3, 9231,
                                                                       9237, 22149, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37954, 0, 3,
                                                                       37759, 22029, 37774, 9249,
                                                                       9267, 22159, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37999, 0, 3,
                                                                       37774, 22039, 37789, 9267,
                                                                       9285, 22189, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38044, 0, 3,
                                                                       37789, 22049, 37804, 9285,
                                                                       9303, 22219, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38089, 0, 3,
                                                                       37804, 22059, 37819, 9303,
                                                                       9321, 22249, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38134, 0, 3,
                                                                       37819, 22069, 37834, 9321,
                                                                       9339, 22279, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38179, 0, 3,
                                                                       37834, 22079, 37849, 9339,
                                                                       9357, 22309, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38224, 0, 3,
                                                                       37849, 22089, 37864, 9357,
                                                                       9375, 22339, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38269, 0, 3,
                                                                       37864, 22099, 37879, 9375,
                                                                       9393, 22369, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38314, 0, 3,
                                                                       37879, 22109, 37894, 9393,
                                                                       9411, 22399, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38359, 0, 3,
                                                                       37894, 22119, 37909, 9411,
                                                                       9429, 22429, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38404, 0, 3,
                                                                       37909, 22129, 37924, 9429,
                                                                       9447, 22459, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38449, 0, 3,
                                                                       37924, 22139, 37939, 9447,
                                                                       9465, 22489, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38494, 0, 3,
                                                                       37954, 22159, 37999, 9501,
                                                                       9537, 22519, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38584, 0, 3,
                                                                       37999, 22189, 38044, 9537,
                                                                       9573, 22579, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38674, 0, 3,
                                                                       38044, 22219, 38089, 9573,
                                                                       9609, 22639, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38764, 0, 3,
                                                                       38089, 22249, 38134, 9609,
                                                                       9645, 22699, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38854, 0, 3,
                                                                       38134, 22279, 38179, 9645,
                                                                       9681, 22759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38944, 0, 3,
                                                                       38179, 22309, 38224, 9681,
                                                                       9717, 22819, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39034, 0, 3,
                                                                       38224, 22339, 38269, 9717,
                                                                       9753, 22879, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39124, 0, 3,
                                                                       38269, 22369, 38314, 9753,
                                                                       9789, 22939, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39214, 0, 3,
                                                                       38314, 22399, 38359, 9789,
                                                                       9825, 22999, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39304, 0, 3,
                                                                       38359, 22429, 38404, 9825,
                                                                       9861, 23059, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39394, 0, 3,
                                                                       38404, 22459, 38449, 9861,
                                                                       9897, 23119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39484, 0, 3,
                                                                       38494, 22519, 38584, 9969,
                                                                       10029, 23179, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39634, 0, 3,
                                                                       38584, 22579, 38674,
                                                                       10029, 10089, 23279,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39784, 0, 3,
                                                                       38674, 22639, 38764,
                                                                       10089, 10149, 23379,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39934, 0, 3,
                                                                       38764, 22699, 38854,
                                                                       10149, 10209, 23479,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40084, 0, 3,
                                                                       38854, 22759, 38944,
                                                                       10209, 10269, 23579,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40234, 0, 3,
                                                                       38944, 22819, 39034,
                                                                       10269, 10329, 23679,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40384, 0, 3,
                                                                       39034, 22879, 39124,
                                                                       10329, 10389, 23779,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40534, 0, 3,
                                                                       39124, 22939, 39214,
                                                                       10389, 10449, 23879,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40684, 0, 3,
                                                                       39214, 22999, 39304,
                                                                       10449, 10509, 23979,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40834, 0, 3,
                                                                       39304, 23059, 39394,
                                                                       10509, 10569, 24079,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40984, 0, 3,
                                                                       39484, 23179, 39634,
                                                                       10689, 10779, 24179,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41209, 0, 3,
                                                                       39634, 23279, 39784,
                                                                       10779, 10869, 24329,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41434, 0, 3,
                                                                       39784, 23379, 39934,
                                                                       10869, 10959, 24479,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41659, 0, 3,
                                                                       39934, 23479, 40084,
                                                                       10959, 11049, 24629,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41884, 0, 3,
                                                                       40084, 23579, 40234,
                                                                       11049, 11139, 24779,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42109, 0, 3,
                                                                       40234, 23679, 40384,
                                                                       11139, 11229, 24929,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42334, 0, 3,
                                                                       40384, 23779, 40534,
                                                                       11229, 11319, 25079,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42559, 0, 3,
                                                                       40534, 23879, 40684,
                                                                       11319, 11409, 25229,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42784, 0, 3,
                                                                       40684, 23979, 40834,
                                                                       11409, 11499, 25379,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43009, 0, 3,
                                                                       40984, 24179, 41209,
                                                                       11679, 11805, 25529,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43324, 0, 3,
                                                                       41209, 24329, 41434,
                                                                       11805, 11931, 25739,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43639, 0, 3,
                                                                       41434, 24479, 41659,
                                                                       11931, 12057, 25949,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43954, 0, 3,
                                                                       41659, 24629, 41884,
                                                                       12057, 12183, 26159,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44269, 0, 3,
                                                                       41884, 24779, 42109,
                                                                       12183, 12309, 26369,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44584, 0, 3,
                                                                       42109, 24929, 42334,
                                                                       12309, 12435, 26579,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44899, 0, 3,
                                                                       42334, 25079, 42559,
                                                                       12435, 12561, 26789,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 45214, 0, 3,
                                                                       42559, 25229, 42784,
                                                                       12561, 12687, 26999,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 45529, 0, 3,
                                                                       43009, 25529, 43324,
                                                                       12939, 13107, 27209,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 45949, 0, 3,
                                                                       43324, 25739, 43639,
                                                                       13107, 13275, 27489,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46369, 0, 3,
                                                                       43639, 25949, 43954,
                                                                       13275, 13443, 27769,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46789, 0, 3,
                                                                       43954, 26159, 44269,
                                                                       13443, 13611, 28049,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 47209, 0, 3,
                                                                       44269, 26369, 44584,
                                                                       13611, 13779, 28329,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 47629, 0, 3,
                                                                       44584, 26579, 44899,
                                                                       13779, 13947, 28609,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 48049, 0, 3,
                                                                       44899, 26789, 45214,
                                                                       13947, 14115, 28889,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 48469, 0, 3,
                                                                       45529, 27209, 45949,
                                                                       14451, 14667, 29169,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 49009, 0, 3,
                                                                       45949, 27489, 46369,
                                                                       14667, 14883, 29529,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 49549, 0, 3,
                                                                       46369, 27769, 46789,
                                                                       14883, 15099, 29889,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 50089, 0, 3,
                                                                       46789, 28049, 47209,
                                                                       15099, 15315, 30249,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 50629, 0, 3,
                                                                       47209, 28329, 47629,
                                                                       15315, 15531, 30609,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 51169, 0, 3,
                                                                       47629, 28609, 48049,
                                                                       15531, 15747, 30969,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 51709, 0, 3,
                                                                       48469, 29169, 49009,
                                                                       16179, 16449, 31329,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 52384, 0, 3,
                                                                       49009, 29529, 49549,
                                                                       16449, 16719, 31779,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 53059, 0, 3,
                                                                       49549, 29889, 50089,
                                                                       16719, 16989, 32229,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 53734, 0, 3,
                                                                       50089, 30249, 50629,
                                                                       16989, 17259, 32679,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 54409, 0, 3,
                                                                       50629, 30609, 51169,
                                                                       17259, 17529, 33129,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 55084, 0, 3,
                                                                       51709, 31329, 52384,
                                                                       18069, 18399, 33579,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 55909, 0, 3,
                                                                       52384, 31779, 53059,
                                                                       18399, 18729, 34129,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 56734, 0, 3,
                                                                       53059, 32229, 53734,
                                                                       18729, 19059, 34679,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 57559, 0, 3,
                                                                       53734, 32679, 54409,
                                                                       19059, 19389, 35229,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 58384, 0, 3,
                                                                       55084, 33579, 55909,
                                                                       20049, 20445, 35779,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 59374, 0, 3,
                                                                       55909, 34129, 56734,
                                                                       20445, 20841, 36439,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 60364, 0, 3,
                                                                       56734, 34679, 57559,
                                                                       20841, 21237, 37099,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61354, 3, 22029,
                                                                       22039, 37789, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61375, 3, 22039,
                                                                       22049, 37804, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61396, 3, 22049,
                                                                       22059, 37819, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61417, 3, 22059,
                                                                       22069, 37834, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61438, 3, 22069,
                                                                       22079, 37849, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61459, 3, 22079,
                                                                       22089, 37864, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61480, 3, 22089,
                                                                       22099, 37879, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61501, 3, 22099,
                                                                       22109, 37894, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61522, 3, 22109,
                                                                       22119, 37909, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61543, 3, 22119,
                                                                       22129, 37924, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61564, 3, 22129,
                                                                       22139, 37939, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61585, 0, 3,
                                                                       61354, 37789, 61375,
                                                                       22159, 22189, 38044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61648, 0, 3,
                                                                       61375, 37804, 61396,
                                                                       22189, 22219, 38089,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61711, 0, 3,
                                                                       61396, 37819, 61417,
                                                                       22219, 22249, 38134,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61774, 0, 3,
                                                                       61417, 37834, 61438,
                                                                       22249, 22279, 38179,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61837, 0, 3,
                                                                       61438, 37849, 61459,
                                                                       22279, 22309, 38224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61900, 0, 3,
                                                                       61459, 37864, 61480,
                                                                       22309, 22339, 38269,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61963, 0, 3,
                                                                       61480, 37879, 61501,
                                                                       22339, 22369, 38314,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62026, 0, 3,
                                                                       61501, 37894, 61522,
                                                                       22369, 22399, 38359,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62089, 0, 3,
                                                                       61522, 37909, 61543,
                                                                       22399, 22429, 38404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62152, 0, 3,
                                                                       61543, 37924, 61564,
                                                                       22429, 22459, 38449,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62215, 0, 3,
                                                                       61585, 38044, 61648,
                                                                       22519, 22579, 38674,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62341, 0, 3,
                                                                       61648, 38089, 61711,
                                                                       22579, 22639, 38764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62467, 0, 3,
                                                                       61711, 38134, 61774,
                                                                       22639, 22699, 38854,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62593, 0, 3,
                                                                       61774, 38179, 61837,
                                                                       22699, 22759, 38944,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62719, 0, 3,
                                                                       61837, 38224, 61900,
                                                                       22759, 22819, 39034,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62845, 0, 3,
                                                                       61900, 38269, 61963,
                                                                       22819, 22879, 39124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62971, 0, 3,
                                                                       61963, 38314, 62026,
                                                                       22879, 22939, 39214,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 63097, 0, 3,
                                                                       62026, 38359, 62089,
                                                                       22939, 22999, 39304,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 63223, 0, 3,
                                                                       62089, 38404, 62152,
                                                                       22999, 23059, 39394,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63349, 0, 3,
                                                                       62215, 38674, 62341,
                                                                       23179, 23279, 39784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63559, 0, 3,
                                                                       62341, 38764, 62467,
                                                                       23279, 23379, 39934,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63769, 0, 3,
                                                                       62467, 38854, 62593,
                                                                       23379, 23479, 40084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63979, 0, 3,
                                                                       62593, 38944, 62719,
                                                                       23479, 23579, 40234,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64189, 0, 3,
                                                                       62719, 39034, 62845,
                                                                       23579, 23679, 40384,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64399, 0, 3,
                                                                       62845, 39124, 62971,
                                                                       23679, 23779, 40534,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64609, 0, 3,
                                                                       62971, 39214, 63097,
                                                                       23779, 23879, 40684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64819, 0, 3,
                                                                       63097, 39304, 63223,
                                                                       23879, 23979, 40834,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65029, 0, 3,
                                                                       63349, 39784, 63559,
                                                                       24179, 24329, 41434,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65344, 0, 3,
                                                                       63559, 39934, 63769,
                                                                       24329, 24479, 41659,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65659, 0, 3,
                                                                       63769, 40084, 63979,
                                                                       24479, 24629, 41884,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65974, 0, 3,
                                                                       63979, 40234, 64189,
                                                                       24629, 24779, 42109,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 66289, 0, 3,
                                                                       64189, 40384, 64399,
                                                                       24779, 24929, 42334,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 66604, 0, 3,
                                                                       64399, 40534, 64609,
                                                                       24929, 25079, 42559,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 66919, 0, 3,
                                                                       64609, 40684, 64819,
                                                                       25079, 25229, 42784,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 67234, 0, 3,
                                                                       65029, 41434, 65344,
                                                                       25529, 25739, 43639,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 67675, 0, 3,
                                                                       65344, 41659, 65659,
                                                                       25739, 25949, 43954,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 68116, 0, 3,
                                                                       65659, 41884, 65974,
                                                                       25949, 26159, 44269,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 68557, 0, 3,
                                                                       65974, 42109, 66289,
                                                                       26159, 26369, 44584,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 68998, 0, 3,
                                                                       66289, 42334, 66604,
                                                                       26369, 26579, 44899,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 69439, 0, 3,
                                                                       66604, 42559, 66919,
                                                                       26579, 26789, 45214,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 69880, 0, 3,
                                                                       67234, 43639, 67675,
                                                                       27209, 27489, 46369,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 70468, 0, 3,
                                                                       67675, 43954, 68116,
                                                                       27489, 27769, 46789,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 71056, 0, 3,
                                                                       68116, 44269, 68557,
                                                                       27769, 28049, 47209,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 71644, 0, 3,
                                                                       68557, 44584, 68998,
                                                                       28049, 28329, 47629,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 72232, 0, 3,
                                                                       68998, 44899, 69439,
                                                                       28329, 28609, 48049,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 72820, 0, 3,
                                                                       69880, 46369, 70468,
                                                                       29169, 29529, 49549,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 73576, 0, 3,
                                                                       70468, 46789, 71056,
                                                                       29529, 29889, 50089,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 74332, 0, 3,
                                                                       71056, 47209, 71644,
                                                                       29889, 30249, 50629,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 75088, 0, 3,
                                                                       71644, 47629, 72232,
                                                                       30249, 30609, 51169,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 75844, 0, 3,
                                                                       72820, 49549, 73576,
                                                                       31329, 31779, 53059,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 76789, 0, 3,
                                                                       73576, 50089, 74332,
                                                                       31779, 32229, 53734,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 77734, 0, 3,
                                                                       74332, 50629, 75088,
                                                                       32229, 32679, 54409,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 78679, 0, 3,
                                                                       75844, 53059, 76789,
                                                                       33579, 34129, 56734,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 79834, 0, 3,
                                                                       76789, 53734, 77734,
                                                                       34129, 34679, 57559,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 80989, 0, 3,
                                                                       78679, 56734, 79834,
                                                                       35779, 36439, 60364,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82375, 3, 37759,
                                                                       37774, 61354, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82403, 3, 37774,
                                                                       37789, 61375, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82431, 3, 37789,
                                                                       37804, 61396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82459, 3, 37804,
                                                                       37819, 61417, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82487, 3, 37819,
                                                                       37834, 61438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82515, 3, 37834,
                                                                       37849, 61459, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82543, 3, 37849,
                                                                       37864, 61480, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82571, 3, 37864,
                                                                       37879, 61501, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82599, 3, 37879,
                                                                       37894, 61522, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82627, 3, 37894,
                                                                       37909, 61543, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82655, 3, 37909,
                                                                       37924, 61564, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 82683, 0, 3,
                                                                       82375, 61354, 82403,
                                                                       37954, 37999, 61585,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 82767, 0, 3,
                                                                       82403, 61375, 82431,
                                                                       37999, 38044, 61648,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 82851, 0, 3,
                                                                       82431, 61396, 82459,
                                                                       38044, 38089, 61711,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 82935, 0, 3,
                                                                       82459, 61417, 82487,
                                                                       38089, 38134, 61774,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83019, 0, 3,
                                                                       82487, 61438, 82515,
                                                                       38134, 38179, 61837,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83103, 0, 3,
                                                                       82515, 61459, 82543,
                                                                       38179, 38224, 61900,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83187, 0, 3,
                                                                       82543, 61480, 82571,
                                                                       38224, 38269, 61963,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83271, 0, 3,
                                                                       82571, 61501, 82599,
                                                                       38269, 38314, 62026,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83355, 0, 3,
                                                                       82599, 61522, 82627,
                                                                       38314, 38359, 62089,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83439, 0, 3,
                                                                       82627, 61543, 82655,
                                                                       38359, 38404, 62152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 83523, 0, 3,
                                                                       82683, 61585, 82767,
                                                                       38494, 38584, 62215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 83691, 0, 3,
                                                                       82767, 61648, 82851,
                                                                       38584, 38674, 62341,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 83859, 0, 3,
                                                                       82851, 61711, 82935,
                                                                       38674, 38764, 62467,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84027, 0, 3,
                                                                       82935, 61774, 83019,
                                                                       38764, 38854, 62593,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84195, 0, 3,
                                                                       83019, 61837, 83103,
                                                                       38854, 38944, 62719,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84363, 0, 3,
                                                                       83103, 61900, 83187,
                                                                       38944, 39034, 62845,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84531, 0, 3,
                                                                       83187, 61963, 83271,
                                                                       39034, 39124, 62971,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84699, 0, 3,
                                                                       83271, 62026, 83355,
                                                                       39124, 39214, 63097,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84867, 0, 3,
                                                                       83355, 62089, 83439,
                                                                       39214, 39304, 63223,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 85035, 0, 3,
                                                                       83523, 62215, 83691,
                                                                       39484, 39634, 63349,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 85315, 0, 3,
                                                                       83691, 62341, 83859,
                                                                       39634, 39784, 63559,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 85595, 0, 3,
                                                                       83859, 62467, 84027,
                                                                       39784, 39934, 63769,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 85875, 0, 3,
                                                                       84027, 62593, 84195,
                                                                       39934, 40084, 63979,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 86155, 0, 3,
                                                                       84195, 62719, 84363,
                                                                       40084, 40234, 64189,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 86435, 0, 3,
                                                                       84363, 62845, 84531,
                                                                       40234, 40384, 64399,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 86715, 0, 3,
                                                                       84531, 62971, 84699,
                                                                       40384, 40534, 64609,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 86995, 0, 3,
                                                                       84699, 63097, 84867,
                                                                       40534, 40684, 64819,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 87275, 0, 3,
                                                                       85035, 63349, 85315,
                                                                       40984, 41209, 65029,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 87695, 0, 3,
                                                                       85315, 63559, 85595,
                                                                       41209, 41434, 65344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 88115, 0, 3,
                                                                       85595, 63769, 85875,
                                                                       41434, 41659, 65659,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 88535, 0, 3,
                                                                       85875, 63979, 86155,
                                                                       41659, 41884, 65974,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 88955, 0, 3,
                                                                       86155, 64189, 86435,
                                                                       41884, 42109, 66289,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 89375, 0, 3,
                                                                       86435, 64399, 86715,
                                                                       42109, 42334, 66604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 89795, 0, 3,
                                                                       86715, 64609, 86995,
                                                                       42334, 42559, 66919,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 90215, 0, 3,
                                                                       87275, 65029, 87695,
                                                                       43009, 43324, 67234,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 90803, 0, 3,
                                                                       87695, 65344, 88115,
                                                                       43324, 43639, 67675,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 91391, 0, 3,
                                                                       88115, 65659, 88535,
                                                                       43639, 43954, 68116,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 91979, 0, 3,
                                                                       88535, 65974, 88955,
                                                                       43954, 44269, 68557,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 92567, 0, 3,
                                                                       88955, 66289, 89375,
                                                                       44269, 44584, 68998,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 93155, 0, 3,
                                                                       89375, 66604, 89795,
                                                                       44584, 44899, 69439,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 93743, 0, 3,
                                                                       90215, 67234, 90803,
                                                                       45529, 45949, 69880,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 94527, 0, 3,
                                                                       90803, 67675, 91391,
                                                                       45949, 46369, 70468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 95311, 0, 3,
                                                                       91391, 68116, 91979,
                                                                       46369, 46789, 71056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 96095, 0, 3,
                                                                       91979, 68557, 92567,
                                                                       46789, 47209, 71644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 96879, 0, 3,
                                                                       92567, 68998, 93155,
                                                                       47209, 47629, 72232,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 97663, 0, 3,
                                                                       93743, 69880, 94527,
                                                                       48469, 49009, 72820,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 98671, 0, 3,
                                                                       94527, 70468, 95311,
                                                                       49009, 49549, 73576,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 99679, 0, 3,
                                                                       95311, 71056, 96095,
                                                                       49549, 50089, 74332,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 100687, 0, 3,
                                                                       96095, 71644, 96879,
                                                                       50089, 50629, 75088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 101695, 0, 3,
                                                                       97663, 72820, 98671,
                                                                       51709, 52384, 75844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 102955, 0, 3,
                                                                       98671, 73576, 99679,
                                                                       52384, 53059, 76789,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 104215, 0, 3,
                                                                       99679, 74332, 100687,
                                                                       53059, 53734, 77734,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 105475, 0, 3,
                                                                       101695, 75844, 102955,
                                                                       55084, 55909, 78679,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 107015, 0, 3,
                                                                       102955, 76789, 104215,
                                                                       55909, 56734, 79834,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 108555, 0, 3,
                                                                       105475, 78679, 107015,
                                                                       58384, 59374, 80989,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 110403, 93743, 784, ncols);

                    simdfunc::contract_primitives(buffer, 111551, 97663, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 113027, 101695, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 114872, 105475, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 117127, 108555, 1848, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 111187, 110403, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 112559, 111551, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 114287, 113027, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 116412, 114872, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 118975, 117127, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 119833, 111187, 112559, 13, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 120925, 112559, 114287, 13, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 122329, 114287, 116412, 13, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 124084, 116412, 118975, 13, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 126229, 119833, 120925, 13, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 128413, 120925, 122329, 13, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 131221, 122329, 124084, 13, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 134731, 126229, 128413, 13, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 138371, 128413, 131221, 13, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 143051, 134731, 138371, 13, nmax);

        simdtrf::transform_i_inner(buffer, 148511, 143051, 15, 13, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 148511, 169, nmax);
    }

    for (size_t m = 0; m < 1521; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
