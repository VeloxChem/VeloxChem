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


#include "SimdThreeCenterElectronRepulsionRecHIL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSML.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"
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
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hil_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hil_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 416913, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2431 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 416913, 319872, 17770, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 19,
                                                             ncols, fj, mu, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 69, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 72, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 75, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 78, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 81, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 84, 0, 3, 7, 8,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 90, 0, 3, 8, 9,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 96, 0, 3, 9, 10,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 102, 0, 3, 10, 11,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 108, 0, 3, 11, 12,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 114, 0, 3, 12, 13,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 120, 0, 3, 13, 14,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 126, 0, 3, 14, 15,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 132, 0, 3, 15, 16,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 138, 0, 3, 16, 17,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 144, 0, 3, 17, 18,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 150, 0, 3, 18, 19,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 156, 0, 3, 19, 20,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 162, 0, 3, 20, 21,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 168, 0, 3, 21, 22,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 174, 0, 3, 22, 23,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 180, 0, 3, 23, 24,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 186, 0, 3, 24, 25,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 27, 30,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 202, 0, 3, 30, 33,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 212, 0, 3, 33, 36,
                                                                       96, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 36, 39,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 232, 0, 3, 39, 42,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 242, 0, 3, 42, 45,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 252, 0, 3, 45, 48,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 262, 0, 3, 48, 51,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 272, 0, 3, 51, 54,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 282, 0, 3, 54, 57,
                                                                       138, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 292, 0, 3, 57, 60,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 302, 0, 3, 60, 63,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 312, 0, 3, 63, 66,
                                                                       156, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 322, 0, 3, 66, 69,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 332, 0, 3, 69, 72,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 342, 0, 3, 72, 75,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 352, 0, 3, 75, 78,
                                                                       180, 186, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 362, 0, 3, 84, 90,
                                                                       192, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 377, 0, 3, 90, 96,
                                                                       202, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 392, 0, 3, 96,
                                                                       102, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 407, 0, 3, 102,
                                                                       108, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 422, 0, 3, 108,
                                                                       114, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 437, 0, 3, 114,
                                                                       120, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 452, 0, 3, 120,
                                                                       126, 252, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 467, 0, 3, 126,
                                                                       132, 262, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 482, 0, 3, 132,
                                                                       138, 272, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 497, 0, 3, 138,
                                                                       144, 282, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 512, 0, 3, 144,
                                                                       150, 292, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 527, 0, 3, 150,
                                                                       156, 302, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 542, 0, 3, 156,
                                                                       162, 312, 322, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 557, 0, 3, 162,
                                                                       168, 322, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 572, 0, 3, 168,
                                                                       174, 332, 342, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 587, 0, 3, 174,
                                                                       180, 342, 352, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 602, 0, 3, 192,
                                                                       202, 362, 377, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 623, 0, 3, 202,
                                                                       212, 377, 392, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 644, 0, 3, 212,
                                                                       222, 392, 407, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 665, 0, 3, 222,
                                                                       232, 407, 422, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 686, 0, 3, 232,
                                                                       242, 422, 437, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 707, 0, 3, 242,
                                                                       252, 437, 452, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 728, 0, 3, 252,
                                                                       262, 452, 467, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 749, 0, 3, 262,
                                                                       272, 467, 482, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 770, 0, 3, 272,
                                                                       282, 482, 497, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 791, 0, 3, 282,
                                                                       292, 497, 512, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 812, 0, 3, 292,
                                                                       302, 512, 527, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 833, 0, 3, 302,
                                                                       312, 527, 542, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 854, 0, 3, 312,
                                                                       322, 542, 557, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 875, 0, 3, 322,
                                                                       332, 557, 572, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 896, 0, 3, 332,
                                                                       342, 572, 587, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 917, 0, 3, 362,
                                                                       377, 602, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 945, 0, 3, 377,
                                                                       392, 623, 644, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 973, 0, 3, 392,
                                                                       407, 644, 665, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1001, 0, 3, 407,
                                                                       422, 665, 686, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1029, 0, 3, 422,
                                                                       437, 686, 707, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1057, 0, 3, 437,
                                                                       452, 707, 728, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1085, 0, 3, 452,
                                                                       467, 728, 749, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1113, 0, 3, 467,
                                                                       482, 749, 770, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1141, 0, 3, 482,
                                                                       497, 770, 791, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1169, 0, 3, 497,
                                                                       512, 791, 812, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1197, 0, 3, 512,
                                                                       527, 812, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1225, 0, 3, 527,
                                                                       542, 833, 854, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1253, 0, 3, 542,
                                                                       557, 854, 875, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1281, 0, 3, 557,
                                                                       572, 875, 896, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1309, 0, 3, 602,
                                                                       623, 917, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1345, 0, 3, 623,
                                                                       644, 945, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1381, 0, 3, 644,
                                                                       665, 973, 1001, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1417, 0, 3, 665,
                                                                       686, 1001, 1029, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1453, 0, 3, 686,
                                                                       707, 1029, 1057, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1489, 0, 3, 707,
                                                                       728, 1057, 1085, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1525, 0, 3, 728,
                                                                       749, 1085, 1113, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1561, 0, 3, 749,
                                                                       770, 1113, 1141, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1597, 0, 3, 770,
                                                                       791, 1141, 1169, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1633, 0, 3, 791,
                                                                       812, 1169, 1197, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1669, 0, 3, 812,
                                                                       833, 1197, 1225, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1705, 0, 3, 833,
                                                                       854, 1225, 1253, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1741, 0, 3, 854,
                                                                       875, 1253, 1281, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1777, 0, 3, 917,
                                                                       945, 1309, 1345, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1822, 0, 3, 945,
                                                                       973, 1345, 1381, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1867, 0, 3, 973,
                                                                       1001, 1381, 1417, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 1001,
                                                                       1029, 1417, 1453, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1957, 0, 3, 1029,
                                                                       1057, 1453, 1489, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2002, 0, 3, 1057,
                                                                       1085, 1489, 1525, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2047, 0, 3, 1085,
                                                                       1113, 1525, 1561, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2092, 0, 3, 1113,
                                                                       1141, 1561, 1597, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2137, 0, 3, 1141,
                                                                       1169, 1597, 1633, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2182, 0, 3, 1169,
                                                                       1197, 1633, 1669, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2227, 0, 3, 1197,
                                                                       1225, 1669, 1705, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2272, 0, 3, 1225,
                                                                       1253, 1705, 1741, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2317, 0, 3, 1309,
                                                                       1345, 1777, 1822, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1345,
                                                                       1381, 1822, 1867, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2427, 0, 3, 1381,
                                                                       1417, 1867, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2482, 0, 3, 1417,
                                                                       1453, 1912, 1957, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2537, 0, 3, 1453,
                                                                       1489, 1957, 2002, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2592, 0, 3, 1489,
                                                                       1525, 2002, 2047, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2647, 0, 3, 1525,
                                                                       1561, 2047, 2092, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2702, 0, 3, 1561,
                                                                       1597, 2092, 2137, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2757, 0, 3, 1597,
                                                                       1633, 2137, 2182, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2812, 0, 3, 1633,
                                                                       1669, 2182, 2227, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2867, 0, 3, 1669,
                                                                       1705, 2227, 2272, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2922, 0, 3, 1777,
                                                                       1822, 2317, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2988, 0, 3, 1822,
                                                                       1867, 2372, 2427, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3054, 0, 3, 1867,
                                                                       1912, 2427, 2482, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3120, 0, 3, 1912,
                                                                       1957, 2482, 2537, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3186, 0, 3, 1957,
                                                                       2002, 2537, 2592, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3252, 0, 3, 2002,
                                                                       2047, 2592, 2647, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3318, 0, 3, 2047,
                                                                       2092, 2647, 2702, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3384, 0, 3, 2092,
                                                                       2137, 2702, 2757, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3450, 0, 3, 2137,
                                                                       2182, 2757, 2812, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3516, 0, 3, 2182,
                                                                       2227, 2812, 2867, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3582, 0, 3, 2317,
                                                                       2372, 2922, 2988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3660, 0, 3, 2372,
                                                                       2427, 2988, 3054, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3738, 0, 3, 2427,
                                                                       2482, 3054, 3120, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3816, 0, 3, 2482,
                                                                       2537, 3120, 3186, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3894, 0, 3, 2537,
                                                                       2592, 3186, 3252, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3972, 0, 3, 2592,
                                                                       2647, 3252, 3318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 4050, 0, 3, 2647,
                                                                       2702, 3318, 3384, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 4128, 0, 3, 2702,
                                                                       2757, 3384, 3450, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 4206, 0, 3, 2757,
                                                                       2812, 3450, 3516, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4284, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4287, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4290, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4293, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4296, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4299, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4302, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4305, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4308, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4311, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4314, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4317, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4320, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4323, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4326, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4329, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4332, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4335, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4338, 3, 9, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4347, 3, 10, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4356, 3, 11, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4365, 3, 12, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4374, 3, 13, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4383, 3, 14, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4392, 3, 15, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4401, 3, 16, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4410, 3, 17, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4419, 3, 18, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4428, 3, 19, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4437, 3, 20, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4446, 3, 21, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4455, 3, 22, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4464, 3, 23, 75,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4473, 3, 24, 78,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4482, 3, 25, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4491, 3, 33, 96,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4509, 3, 36, 102,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4527, 3, 39, 108,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4545, 3, 42, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4563, 3, 45, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4581, 3, 48, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4599, 3, 51, 132,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4617, 3, 54, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4635, 3, 57, 144,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4653, 3, 60, 150,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4671, 3, 63, 156,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4689, 3, 66, 162,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4707, 3, 69, 168,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4725, 3, 72, 174,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4743, 3, 75, 180,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4761, 3, 78, 186,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4779, 3, 96, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4809, 3, 102, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4839, 3, 108, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4869, 3, 114, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4899, 3, 120, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4929, 3, 126, 262,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4959, 3, 132, 272,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4989, 3, 138, 282,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5019, 3, 144, 292,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5049, 3, 150, 302,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5079, 3, 156, 312,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5109, 3, 162, 322,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5139, 3, 168, 332,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5169, 3, 174, 342,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5199, 3, 180, 352,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5229, 3, 212, 392,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5274, 3, 222, 407,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5319, 3, 232, 422,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5364, 3, 242, 437,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5409, 3, 252, 452,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5454, 3, 262, 467,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5499, 3, 272, 482,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5544, 3, 282, 497,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5589, 3, 292, 512,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5634, 3, 302, 527,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5679, 3, 312, 542,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5724, 3, 322, 557,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5769, 3, 332, 572,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5814, 3, 342, 587,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5859, 3, 392, 644,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5922, 3, 407, 665,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5985, 3, 422, 686,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6048, 3, 437, 707,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6111, 3, 452, 728,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6174, 3, 467, 749,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6237, 3, 482, 770,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6300, 3, 497, 791,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6363, 3, 512, 812,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6426, 3, 527, 833,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6489, 3, 542, 854,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6552, 3, 557, 875,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6615, 3, 572, 896,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6678, 3, 644, 973,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6762, 3, 665,
                                                                       1001, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6846, 3, 686,
                                                                       1029, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6930, 3, 707,
                                                                       1057, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7014, 3, 728,
                                                                       1085, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7098, 3, 749,
                                                                       1113, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7182, 3, 770,
                                                                       1141, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7266, 3, 791,
                                                                       1169, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7350, 3, 812,
                                                                       1197, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7434, 3, 833,
                                                                       1225, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7518, 3, 854,
                                                                       1253, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7602, 3, 875,
                                                                       1281, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7686, 3, 973,
                                                                       1381, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7794, 3, 1001,
                                                                       1417, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7902, 3, 1029,
                                                                       1453, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8010, 3, 1057,
                                                                       1489, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8118, 3, 1085,
                                                                       1525, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8226, 3, 1113,
                                                                       1561, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8334, 3, 1141,
                                                                       1597, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8442, 3, 1169,
                                                                       1633, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8550, 3, 1197,
                                                                       1669, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8658, 3, 1225,
                                                                       1705, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8766, 3, 1253,
                                                                       1741, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8874, 3, 1381,
                                                                       1867, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9009, 3, 1417,
                                                                       1912, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9144, 3, 1453,
                                                                       1957, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9279, 3, 1489,
                                                                       2002, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9414, 3, 1525,
                                                                       2047, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9549, 3, 1561,
                                                                       2092, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9684, 3, 1597,
                                                                       2137, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9819, 3, 1633,
                                                                       2182, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9954, 3, 1669,
                                                                       2227, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10089, 3, 1705,
                                                                       2272, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10224, 3, 1867,
                                                                       2427, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10389, 3, 1912,
                                                                       2482, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10554, 3, 1957,
                                                                       2537, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10719, 3, 2002,
                                                                       2592, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10884, 3, 2047,
                                                                       2647, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11049, 3, 2092,
                                                                       2702, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11214, 3, 2137,
                                                                       2757, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11379, 3, 2182,
                                                                       2812, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11544, 3, 2227,
                                                                       2867, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11709, 3, 2427,
                                                                       3054, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11907, 3, 2482,
                                                                       3120, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12105, 3, 2537,
                                                                       3186, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12303, 3, 2592,
                                                                       3252, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12501, 3, 2647,
                                                                       3318, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12699, 3, 2702,
                                                                       3384, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12897, 3, 2757,
                                                                       3450, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 13095, 3, 2812,
                                                                       3516, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13293, 3, 3054,
                                                                       3738, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13527, 3, 3120,
                                                                       3816, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13761, 3, 3186,
                                                                       3894, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13995, 3, 3252,
                                                                       3972, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 14229, 3, 3318,
                                                                       4050, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 14463, 3, 3384,
                                                                       4128, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 14697, 3, 3450,
                                                                       4206, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14931, 3, 7, 8,
                                                                       4284, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14937, 3, 8, 9,
                                                                       4287, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14943, 3, 9, 10,
                                                                       4290, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14949, 3, 10, 11,
                                                                       4293, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14955, 3, 11, 12,
                                                                       4296, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14961, 3, 12, 13,
                                                                       4299, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14967, 3, 13, 14,
                                                                       4302, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14973, 3, 14, 15,
                                                                       4305, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14979, 3, 15, 16,
                                                                       4308, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14985, 3, 16, 17,
                                                                       4311, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14991, 3, 17, 18,
                                                                       4314, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14997, 3, 18, 19,
                                                                       4317, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15003, 3, 19, 20,
                                                                       4320, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15009, 3, 20, 21,
                                                                       4323, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15015, 3, 21, 22,
                                                                       4326, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15021, 3, 22, 23,
                                                                       4329, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15027, 3, 23, 24,
                                                                       4332, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15033, 3, 24, 25,
                                                                       4335, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15039, 0, 3,
                                                                       14931, 4284, 14937, 4338,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15057, 0, 3,
                                                                       14937, 4287, 14943, 4347,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15075, 0, 3,
                                                                       14943, 4290, 14949, 4356,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15093, 0, 3,
                                                                       14949, 4293, 14955, 4365,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15111, 0, 3,
                                                                       14955, 4296, 14961, 4374,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15129, 0, 3,
                                                                       14961, 4299, 14967, 4383,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15147, 0, 3,
                                                                       14967, 4302, 14973, 4392,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15165, 0, 3,
                                                                       14973, 4305, 14979, 4401,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15183, 0, 3,
                                                                       14979, 4308, 14985, 4410,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15201, 0, 3,
                                                                       14985, 4311, 14991, 4419,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15219, 0, 3,
                                                                       14991, 4314, 14997, 4428,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15237, 0, 3,
                                                                       14997, 4317, 15003, 4437,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15255, 0, 3,
                                                                       15003, 4320, 15009, 4446,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15273, 0, 3,
                                                                       15009, 4323, 15015, 4455,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15291, 0, 3,
                                                                       15015, 4326, 15021, 4464,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15309, 0, 3,
                                                                       15021, 4329, 15027, 4473,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15327, 0, 3,
                                                                       15027, 4332, 15033, 4482,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15345, 0, 3,
                                                                       15039, 4338, 15057, 84,
                                                                       90, 4491, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15381, 0, 3,
                                                                       15057, 4347, 15075, 90,
                                                                       96, 4509, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15417, 0, 3,
                                                                       15075, 4356, 15093, 96,
                                                                       102, 4527, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15453, 0, 3,
                                                                       15093, 4365, 15111, 102,
                                                                       108, 4545, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15489, 0, 3,
                                                                       15111, 4374, 15129, 108,
                                                                       114, 4563, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15525, 0, 3,
                                                                       15129, 4383, 15147, 114,
                                                                       120, 4581, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15561, 0, 3,
                                                                       15147, 4392, 15165, 120,
                                                                       126, 4599, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15597, 0, 3,
                                                                       15165, 4401, 15183, 126,
                                                                       132, 4617, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15633, 0, 3,
                                                                       15183, 4410, 15201, 132,
                                                                       138, 4635, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15669, 0, 3,
                                                                       15201, 4419, 15219, 138,
                                                                       144, 4653, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15705, 0, 3,
                                                                       15219, 4428, 15237, 144,
                                                                       150, 4671, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15741, 0, 3,
                                                                       15237, 4437, 15255, 150,
                                                                       156, 4689, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15777, 0, 3,
                                                                       15255, 4446, 15273, 156,
                                                                       162, 4707, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15813, 0, 3,
                                                                       15273, 4455, 15291, 162,
                                                                       168, 4725, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15849, 0, 3,
                                                                       15291, 4464, 15309, 168,
                                                                       174, 4743, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15885, 0, 3,
                                                                       15309, 4473, 15327, 174,
                                                                       180, 4761, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15921, 0, 3,
                                                                       15345, 4491, 15381, 192,
                                                                       202, 4779, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15981, 0, 3,
                                                                       15381, 4509, 15417, 202,
                                                                       212, 4809, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16041, 0, 3,
                                                                       15417, 4527, 15453, 212,
                                                                       222, 4839, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16101, 0, 3,
                                                                       15453, 4545, 15489, 222,
                                                                       232, 4869, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16161, 0, 3,
                                                                       15489, 4563, 15525, 232,
                                                                       242, 4899, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16221, 0, 3,
                                                                       15525, 4581, 15561, 242,
                                                                       252, 4929, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16281, 0, 3,
                                                                       15561, 4599, 15597, 252,
                                                                       262, 4959, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16341, 0, 3,
                                                                       15597, 4617, 15633, 262,
                                                                       272, 4989, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16401, 0, 3,
                                                                       15633, 4635, 15669, 272,
                                                                       282, 5019, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16461, 0, 3,
                                                                       15669, 4653, 15705, 282,
                                                                       292, 5049, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16521, 0, 3,
                                                                       15705, 4671, 15741, 292,
                                                                       302, 5079, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16581, 0, 3,
                                                                       15741, 4689, 15777, 302,
                                                                       312, 5109, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16641, 0, 3,
                                                                       15777, 4707, 15813, 312,
                                                                       322, 5139, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16701, 0, 3,
                                                                       15813, 4725, 15849, 322,
                                                                       332, 5169, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16761, 0, 3,
                                                                       15849, 4743, 15885, 332,
                                                                       342, 5199, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16821, 0, 3,
                                                                       15921, 4779, 15981, 362,
                                                                       377, 5229, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16911, 0, 3,
                                                                       15981, 4809, 16041, 377,
                                                                       392, 5274, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17001, 0, 3,
                                                                       16041, 4839, 16101, 392,
                                                                       407, 5319, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17091, 0, 3,
                                                                       16101, 4869, 16161, 407,
                                                                       422, 5364, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17181, 0, 3,
                                                                       16161, 4899, 16221, 422,
                                                                       437, 5409, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17271, 0, 3,
                                                                       16221, 4929, 16281, 437,
                                                                       452, 5454, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17361, 0, 3,
                                                                       16281, 4959, 16341, 452,
                                                                       467, 5499, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17451, 0, 3,
                                                                       16341, 4989, 16401, 467,
                                                                       482, 5544, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17541, 0, 3,
                                                                       16401, 5019, 16461, 482,
                                                                       497, 5589, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17631, 0, 3,
                                                                       16461, 5049, 16521, 497,
                                                                       512, 5634, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17721, 0, 3,
                                                                       16521, 5079, 16581, 512,
                                                                       527, 5679, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17811, 0, 3,
                                                                       16581, 5109, 16641, 527,
                                                                       542, 5724, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17901, 0, 3,
                                                                       16641, 5139, 16701, 542,
                                                                       557, 5769, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17991, 0, 3,
                                                                       16701, 5169, 16761, 557,
                                                                       572, 5814, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18081, 0, 3,
                                                                       16821, 5229, 16911, 602,
                                                                       623, 5859, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18207, 0, 3,
                                                                       16911, 5274, 17001, 623,
                                                                       644, 5922, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18333, 0, 3,
                                                                       17001, 5319, 17091, 644,
                                                                       665, 5985, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18459, 0, 3,
                                                                       17091, 5364, 17181, 665,
                                                                       686, 6048, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18585, 0, 3,
                                                                       17181, 5409, 17271, 686,
                                                                       707, 6111, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18711, 0, 3,
                                                                       17271, 5454, 17361, 707,
                                                                       728, 6174, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18837, 0, 3,
                                                                       17361, 5499, 17451, 728,
                                                                       749, 6237, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18963, 0, 3,
                                                                       17451, 5544, 17541, 749,
                                                                       770, 6300, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19089, 0, 3,
                                                                       17541, 5589, 17631, 770,
                                                                       791, 6363, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19215, 0, 3,
                                                                       17631, 5634, 17721, 791,
                                                                       812, 6426, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19341, 0, 3,
                                                                       17721, 5679, 17811, 812,
                                                                       833, 6489, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19467, 0, 3,
                                                                       17811, 5724, 17901, 833,
                                                                       854, 6552, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19593, 0, 3,
                                                                       17901, 5769, 17991, 854,
                                                                       875, 6615, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19719, 0, 3,
                                                                       18081, 5859, 18207, 917,
                                                                       945, 6678, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19887, 0, 3,
                                                                       18207, 5922, 18333, 945,
                                                                       973, 6762, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20055, 0, 3,
                                                                       18333, 5985, 18459, 973,
                                                                       1001, 6846, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20223, 0, 3,
                                                                       18459, 6048, 18585, 1001,
                                                                       1029, 6930, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20391, 0, 3,
                                                                       18585, 6111, 18711, 1029,
                                                                       1057, 7014, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20559, 0, 3,
                                                                       18711, 6174, 18837, 1057,
                                                                       1085, 7098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20727, 0, 3,
                                                                       18837, 6237, 18963, 1085,
                                                                       1113, 7182, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20895, 0, 3,
                                                                       18963, 6300, 19089, 1113,
                                                                       1141, 7266, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21063, 0, 3,
                                                                       19089, 6363, 19215, 1141,
                                                                       1169, 7350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21231, 0, 3,
                                                                       19215, 6426, 19341, 1169,
                                                                       1197, 7434, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21399, 0, 3,
                                                                       19341, 6489, 19467, 1197,
                                                                       1225, 7518, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21567, 0, 3,
                                                                       19467, 6552, 19593, 1225,
                                                                       1253, 7602, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21735, 0, 3,
                                                                       19719, 6678, 19887, 1309,
                                                                       1345, 7686, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21951, 0, 3,
                                                                       19887, 6762, 20055, 1345,
                                                                       1381, 7794, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22167, 0, 3,
                                                                       20055, 6846, 20223, 1381,
                                                                       1417, 7902, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22383, 0, 3,
                                                                       20223, 6930, 20391, 1417,
                                                                       1453, 8010, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22599, 0, 3,
                                                                       20391, 7014, 20559, 1453,
                                                                       1489, 8118, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22815, 0, 3,
                                                                       20559, 7098, 20727, 1489,
                                                                       1525, 8226, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23031, 0, 3,
                                                                       20727, 7182, 20895, 1525,
                                                                       1561, 8334, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23247, 0, 3,
                                                                       20895, 7266, 21063, 1561,
                                                                       1597, 8442, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23463, 0, 3,
                                                                       21063, 7350, 21231, 1597,
                                                                       1633, 8550, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23679, 0, 3,
                                                                       21231, 7434, 21399, 1633,
                                                                       1669, 8658, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23895, 0, 3,
                                                                       21399, 7518, 21567, 1669,
                                                                       1705, 8766, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24111, 0, 3,
                                                                       21735, 7686, 21951, 1777,
                                                                       1822, 8874, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24381, 0, 3,
                                                                       21951, 7794, 22167, 1822,
                                                                       1867, 9009, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24651, 0, 3,
                                                                       22167, 7902, 22383, 1867,
                                                                       1912, 9144, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24921, 0, 3,
                                                                       22383, 8010, 22599, 1912,
                                                                       1957, 9279, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 25191, 0, 3,
                                                                       22599, 8118, 22815, 1957,
                                                                       2002, 9414, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 25461, 0, 3,
                                                                       22815, 8226, 23031, 2002,
                                                                       2047, 9549, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 25731, 0, 3,
                                                                       23031, 8334, 23247, 2047,
                                                                       2092, 9684, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26001, 0, 3,
                                                                       23247, 8442, 23463, 2092,
                                                                       2137, 9819, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26271, 0, 3,
                                                                       23463, 8550, 23679, 2137,
                                                                       2182, 9954, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26541, 0, 3,
                                                                       23679, 8658, 23895, 2182,
                                                                       2227, 10089, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 26811, 0, 3,
                                                                       24111, 8874, 24381, 2317,
                                                                       2372, 10224, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 27141, 0, 3,
                                                                       24381, 9009, 24651, 2372,
                                                                       2427, 10389, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 27471, 0, 3,
                                                                       24651, 9144, 24921, 2427,
                                                                       2482, 10554, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 27801, 0, 3,
                                                                       24921, 9279, 25191, 2482,
                                                                       2537, 10719, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 28131, 0, 3,
                                                                       25191, 9414, 25461, 2537,
                                                                       2592, 10884, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 28461, 0, 3,
                                                                       25461, 9549, 25731, 2592,
                                                                       2647, 11049, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 28791, 0, 3,
                                                                       25731, 9684, 26001, 2647,
                                                                       2702, 11214, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 29121, 0, 3,
                                                                       26001, 9819, 26271, 2702,
                                                                       2757, 11379, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 29451, 0, 3,
                                                                       26271, 9954, 26541, 2757,
                                                                       2812, 11544, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 29781, 0, 3,
                                                                       26811, 10224, 27141, 2922,
                                                                       2988, 11709, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 30177, 0, 3,
                                                                       27141, 10389, 27471, 2988,
                                                                       3054, 11907, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 30573, 0, 3,
                                                                       27471, 10554, 27801, 3054,
                                                                       3120, 12105, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 30969, 0, 3,
                                                                       27801, 10719, 28131, 3120,
                                                                       3186, 12303, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 31365, 0, 3,
                                                                       28131, 10884, 28461, 3186,
                                                                       3252, 12501, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 31761, 0, 3,
                                                                       28461, 11049, 28791, 3252,
                                                                       3318, 12699, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 32157, 0, 3,
                                                                       28791, 11214, 29121, 3318,
                                                                       3384, 12897, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 32553, 0, 3,
                                                                       29121, 11379, 29451, 3384,
                                                                       3450, 13095, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 32949, 0, 3,
                                                                       29781, 11709, 30177, 3582,
                                                                       3660, 13293, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 33417, 0, 3,
                                                                       30177, 11907, 30573, 3660,
                                                                       3738, 13527, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 33885, 0, 3,
                                                                       30573, 12105, 30969, 3738,
                                                                       3816, 13761, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 34353, 0, 3,
                                                                       30969, 12303, 31365, 3816,
                                                                       3894, 13995, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 34821, 0, 3,
                                                                       31365, 12501, 31761, 3894,
                                                                       3972, 14229, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 35289, 0, 3,
                                                                       31761, 12699, 32157, 3972,
                                                                       4050, 14463, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 35757, 0, 3,
                                                                       32157, 12897, 32553, 4050,
                                                                       4128, 14697, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36225, 3, 4284,
                                                                       4287, 14943, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36235, 3, 4287,
                                                                       4290, 14949, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36245, 3, 4290,
                                                                       4293, 14955, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36255, 3, 4293,
                                                                       4296, 14961, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36265, 3, 4296,
                                                                       4299, 14967, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36275, 3, 4299,
                                                                       4302, 14973, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36285, 3, 4302,
                                                                       4305, 14979, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36295, 3, 4305,
                                                                       4308, 14985, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36305, 3, 4308,
                                                                       4311, 14991, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36315, 3, 4311,
                                                                       4314, 14997, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36325, 3, 4314,
                                                                       4317, 15003, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36335, 3, 4317,
                                                                       4320, 15009, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36345, 3, 4320,
                                                                       4323, 15015, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36355, 3, 4323,
                                                                       4326, 15021, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36365, 3, 4326,
                                                                       4329, 15027, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36375, 3, 4329,
                                                                       4332, 15033, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36385, 0, 3,
                                                                       36225, 14943, 36235,
                                                                       15075, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36415, 0, 3,
                                                                       36235, 14949, 36245,
                                                                       15093, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36445, 0, 3,
                                                                       36245, 14955, 36255,
                                                                       15111, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36475, 0, 3,
                                                                       36255, 14961, 36265,
                                                                       15129, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36505, 0, 3,
                                                                       36265, 14967, 36275,
                                                                       15147, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36535, 0, 3,
                                                                       36275, 14973, 36285,
                                                                       15165, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36565, 0, 3,
                                                                       36285, 14979, 36295,
                                                                       15183, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36595, 0, 3,
                                                                       36295, 14985, 36305,
                                                                       15201, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36625, 0, 3,
                                                                       36305, 14991, 36315,
                                                                       15219, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36655, 0, 3,
                                                                       36315, 14997, 36325,
                                                                       15237, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36685, 0, 3,
                                                                       36325, 15003, 36335,
                                                                       15255, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36715, 0, 3,
                                                                       36335, 15009, 36345,
                                                                       15273, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36745, 0, 3,
                                                                       36345, 15015, 36355,
                                                                       15291, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36775, 0, 3,
                                                                       36355, 15021, 36365,
                                                                       15309, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36805, 0, 3,
                                                                       36365, 15027, 36375,
                                                                       15327, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36835, 0, 3,
                                                                       36385, 15075, 36415, 4491,
                                                                       4509, 15417, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36895, 0, 3,
                                                                       36415, 15093, 36445, 4509,
                                                                       4527, 15453, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36955, 0, 3,
                                                                       36445, 15111, 36475, 4527,
                                                                       4545, 15489, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37015, 0, 3,
                                                                       36475, 15129, 36505, 4545,
                                                                       4563, 15525, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37075, 0, 3,
                                                                       36505, 15147, 36535, 4563,
                                                                       4581, 15561, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37135, 0, 3,
                                                                       36535, 15165, 36565, 4581,
                                                                       4599, 15597, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37195, 0, 3,
                                                                       36565, 15183, 36595, 4599,
                                                                       4617, 15633, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37255, 0, 3,
                                                                       36595, 15201, 36625, 4617,
                                                                       4635, 15669, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37315, 0, 3,
                                                                       36625, 15219, 36655, 4635,
                                                                       4653, 15705, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37375, 0, 3,
                                                                       36655, 15237, 36685, 4653,
                                                                       4671, 15741, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37435, 0, 3,
                                                                       36685, 15255, 36715, 4671,
                                                                       4689, 15777, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37495, 0, 3,
                                                                       36715, 15273, 36745, 4689,
                                                                       4707, 15813, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37555, 0, 3,
                                                                       36745, 15291, 36775, 4707,
                                                                       4725, 15849, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37615, 0, 3,
                                                                       36775, 15309, 36805, 4725,
                                                                       4743, 15885, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 37675, 0, 3,
                                                                       36835, 15417, 36895, 4779,
                                                                       4809, 16041, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 37775, 0, 3,
                                                                       36895, 15453, 36955, 4809,
                                                                       4839, 16101, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 37875, 0, 3,
                                                                       36955, 15489, 37015, 4839,
                                                                       4869, 16161, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 37975, 0, 3,
                                                                       37015, 15525, 37075, 4869,
                                                                       4899, 16221, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38075, 0, 3,
                                                                       37075, 15561, 37135, 4899,
                                                                       4929, 16281, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38175, 0, 3,
                                                                       37135, 15597, 37195, 4929,
                                                                       4959, 16341, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38275, 0, 3,
                                                                       37195, 15633, 37255, 4959,
                                                                       4989, 16401, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38375, 0, 3,
                                                                       37255, 15669, 37315, 4989,
                                                                       5019, 16461, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38475, 0, 3,
                                                                       37315, 15705, 37375, 5019,
                                                                       5049, 16521, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38575, 0, 3,
                                                                       37375, 15741, 37435, 5049,
                                                                       5079, 16581, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38675, 0, 3,
                                                                       37435, 15777, 37495, 5079,
                                                                       5109, 16641, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38775, 0, 3,
                                                                       37495, 15813, 37555, 5109,
                                                                       5139, 16701, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38875, 0, 3,
                                                                       37555, 15849, 37615, 5139,
                                                                       5169, 16761, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38975, 0, 3,
                                                                       37675, 16041, 37775, 5229,
                                                                       5274, 17001, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39125, 0, 3,
                                                                       37775, 16101, 37875, 5274,
                                                                       5319, 17091, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39275, 0, 3,
                                                                       37875, 16161, 37975, 5319,
                                                                       5364, 17181, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39425, 0, 3,
                                                                       37975, 16221, 38075, 5364,
                                                                       5409, 17271, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39575, 0, 3,
                                                                       38075, 16281, 38175, 5409,
                                                                       5454, 17361, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39725, 0, 3,
                                                                       38175, 16341, 38275, 5454,
                                                                       5499, 17451, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39875, 0, 3,
                                                                       38275, 16401, 38375, 5499,
                                                                       5544, 17541, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40025, 0, 3,
                                                                       38375, 16461, 38475, 5544,
                                                                       5589, 17631, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40175, 0, 3,
                                                                       38475, 16521, 38575, 5589,
                                                                       5634, 17721, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40325, 0, 3,
                                                                       38575, 16581, 38675, 5634,
                                                                       5679, 17811, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40475, 0, 3,
                                                                       38675, 16641, 38775, 5679,
                                                                       5724, 17901, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40625, 0, 3,
                                                                       38775, 16701, 38875, 5724,
                                                                       5769, 17991, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40775, 0, 3,
                                                                       38975, 17001, 39125, 5859,
                                                                       5922, 18333, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40985, 0, 3,
                                                                       39125, 17091, 39275, 5922,
                                                                       5985, 18459, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41195, 0, 3,
                                                                       39275, 17181, 39425, 5985,
                                                                       6048, 18585, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41405, 0, 3,
                                                                       39425, 17271, 39575, 6048,
                                                                       6111, 18711, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41615, 0, 3,
                                                                       39575, 17361, 39725, 6111,
                                                                       6174, 18837, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41825, 0, 3,
                                                                       39725, 17451, 39875, 6174,
                                                                       6237, 18963, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42035, 0, 3,
                                                                       39875, 17541, 40025, 6237,
                                                                       6300, 19089, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42245, 0, 3,
                                                                       40025, 17631, 40175, 6300,
                                                                       6363, 19215, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42455, 0, 3,
                                                                       40175, 17721, 40325, 6363,
                                                                       6426, 19341, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42665, 0, 3,
                                                                       40325, 17811, 40475, 6426,
                                                                       6489, 19467, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42875, 0, 3,
                                                                       40475, 17901, 40625, 6489,
                                                                       6552, 19593, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43085, 0, 3,
                                                                       40775, 18333, 40985, 6678,
                                                                       6762, 20055, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43365, 0, 3,
                                                                       40985, 18459, 41195, 6762,
                                                                       6846, 20223, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43645, 0, 3,
                                                                       41195, 18585, 41405, 6846,
                                                                       6930, 20391, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43925, 0, 3,
                                                                       41405, 18711, 41615, 6930,
                                                                       7014, 20559, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44205, 0, 3,
                                                                       41615, 18837, 41825, 7014,
                                                                       7098, 20727, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44485, 0, 3,
                                                                       41825, 18963, 42035, 7098,
                                                                       7182, 20895, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44765, 0, 3,
                                                                       42035, 19089, 42245, 7182,
                                                                       7266, 21063, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 45045, 0, 3,
                                                                       42245, 19215, 42455, 7266,
                                                                       7350, 21231, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 45325, 0, 3,
                                                                       42455, 19341, 42665, 7350,
                                                                       7434, 21399, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 45605, 0, 3,
                                                                       42665, 19467, 42875, 7434,
                                                                       7518, 21567, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 45885, 0, 3,
                                                                       43085, 20055, 43365, 7686,
                                                                       7794, 22167, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46245, 0, 3,
                                                                       43365, 20223, 43645, 7794,
                                                                       7902, 22383, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46605, 0, 3,
                                                                       43645, 20391, 43925, 7902,
                                                                       8010, 22599, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46965, 0, 3,
                                                                       43925, 20559, 44205, 8010,
                                                                       8118, 22815, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 47325, 0, 3,
                                                                       44205, 20727, 44485, 8118,
                                                                       8226, 23031, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 47685, 0, 3,
                                                                       44485, 20895, 44765, 8226,
                                                                       8334, 23247, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48045, 0, 3,
                                                                       44765, 21063, 45045, 8334,
                                                                       8442, 23463, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48405, 0, 3,
                                                                       45045, 21231, 45325, 8442,
                                                                       8550, 23679, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48765, 0, 3,
                                                                       45325, 21399, 45605, 8550,
                                                                       8658, 23895, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 49125, 0, 3,
                                                                       45885, 22167, 46245, 8874,
                                                                       9009, 24651, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 49575, 0, 3,
                                                                       46245, 22383, 46605, 9009,
                                                                       9144, 24921, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 50025, 0, 3,
                                                                       46605, 22599, 46965, 9144,
                                                                       9279, 25191, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 50475, 0, 3,
                                                                       46965, 22815, 47325, 9279,
                                                                       9414, 25461, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 50925, 0, 3,
                                                                       47325, 23031, 47685, 9414,
                                                                       9549, 25731, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 51375, 0, 3,
                                                                       47685, 23247, 48045, 9549,
                                                                       9684, 26001, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 51825, 0, 3,
                                                                       48045, 23463, 48405, 9684,
                                                                       9819, 26271, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 52275, 0, 3,
                                                                       48405, 23679, 48765, 9819,
                                                                       9954, 26541, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 52725, 0, 3,
                                                                       49125, 24651, 49575,
                                                                       10224, 10389, 27471,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 53275, 0, 3,
                                                                       49575, 24921, 50025,
                                                                       10389, 10554, 27801,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 53825, 0, 3,
                                                                       50025, 25191, 50475,
                                                                       10554, 10719, 28131,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 54375, 0, 3,
                                                                       50475, 25461, 50925,
                                                                       10719, 10884, 28461,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 54925, 0, 3,
                                                                       50925, 25731, 51375,
                                                                       10884, 11049, 28791,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 55475, 0, 3,
                                                                       51375, 26001, 51825,
                                                                       11049, 11214, 29121,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 56025, 0, 3,
                                                                       51825, 26271, 52275,
                                                                       11214, 11379, 29451,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 56575, 0, 3,
                                                                       52725, 27471, 53275,
                                                                       11709, 11907, 30573,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 57235, 0, 3,
                                                                       53275, 27801, 53825,
                                                                       11907, 12105, 30969,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 57895, 0, 3,
                                                                       53825, 28131, 54375,
                                                                       12105, 12303, 31365,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 58555, 0, 3,
                                                                       54375, 28461, 54925,
                                                                       12303, 12501, 31761,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 59215, 0, 3,
                                                                       54925, 28791, 55475,
                                                                       12501, 12699, 32157,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 59875, 0, 3,
                                                                       55475, 29121, 56025,
                                                                       12699, 12897, 32553,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 60535, 0, 3,
                                                                       56575, 30573, 57235,
                                                                       13293, 13527, 33885,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 61315, 0, 3,
                                                                       57235, 30969, 57895,
                                                                       13527, 13761, 34353,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 62095, 0, 3,
                                                                       57895, 31365, 58555,
                                                                       13761, 13995, 34821,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 62875, 0, 3,
                                                                       58555, 31761, 59215,
                                                                       13995, 14229, 35289,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 63655, 0, 3,
                                                                       59215, 32157, 59875,
                                                                       14229, 14463, 35757,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64435, 3, 14931,
                                                                       14937, 36225, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64450, 3, 14937,
                                                                       14943, 36235, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64465, 3, 14943,
                                                                       14949, 36245, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64480, 3, 14949,
                                                                       14955, 36255, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64495, 3, 14955,
                                                                       14961, 36265, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64510, 3, 14961,
                                                                       14967, 36275, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64525, 3, 14967,
                                                                       14973, 36285, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64540, 3, 14973,
                                                                       14979, 36295, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64555, 3, 14979,
                                                                       14985, 36305, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64570, 3, 14985,
                                                                       14991, 36315, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64585, 3, 14991,
                                                                       14997, 36325, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64600, 3, 14997,
                                                                       15003, 36335, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64615, 3, 15003,
                                                                       15009, 36345, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64630, 3, 15009,
                                                                       15015, 36355, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64645, 3, 15015,
                                                                       15021, 36365, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64660, 3, 15021,
                                                                       15027, 36375, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64675, 0, 3,
                                                                       64435, 36225, 64450,
                                                                       15039, 15057, 36385,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64720, 0, 3,
                                                                       64450, 36235, 64465,
                                                                       15057, 15075, 36415,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64765, 0, 3,
                                                                       64465, 36245, 64480,
                                                                       15075, 15093, 36445,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64810, 0, 3,
                                                                       64480, 36255, 64495,
                                                                       15093, 15111, 36475,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64855, 0, 3,
                                                                       64495, 36265, 64510,
                                                                       15111, 15129, 36505,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64900, 0, 3,
                                                                       64510, 36275, 64525,
                                                                       15129, 15147, 36535,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64945, 0, 3,
                                                                       64525, 36285, 64540,
                                                                       15147, 15165, 36565,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64990, 0, 3,
                                                                       64540, 36295, 64555,
                                                                       15165, 15183, 36595,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65035, 0, 3,
                                                                       64555, 36305, 64570,
                                                                       15183, 15201, 36625,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65080, 0, 3,
                                                                       64570, 36315, 64585,
                                                                       15201, 15219, 36655,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65125, 0, 3,
                                                                       64585, 36325, 64600,
                                                                       15219, 15237, 36685,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65170, 0, 3,
                                                                       64600, 36335, 64615,
                                                                       15237, 15255, 36715,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65215, 0, 3,
                                                                       64615, 36345, 64630,
                                                                       15255, 15273, 36745,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65260, 0, 3,
                                                                       64630, 36355, 64645,
                                                                       15273, 15291, 36775,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65305, 0, 3,
                                                                       64645, 36365, 64660,
                                                                       15291, 15309, 36805,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65350, 0, 3,
                                                                       64675, 36385, 64720,
                                                                       15345, 15381, 36835,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65440, 0, 3,
                                                                       64720, 36415, 64765,
                                                                       15381, 15417, 36895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65530, 0, 3,
                                                                       64765, 36445, 64810,
                                                                       15417, 15453, 36955,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65620, 0, 3,
                                                                       64810, 36475, 64855,
                                                                       15453, 15489, 37015,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65710, 0, 3,
                                                                       64855, 36505, 64900,
                                                                       15489, 15525, 37075,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65800, 0, 3,
                                                                       64900, 36535, 64945,
                                                                       15525, 15561, 37135,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65890, 0, 3,
                                                                       64945, 36565, 64990,
                                                                       15561, 15597, 37195,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65980, 0, 3,
                                                                       64990, 36595, 65035,
                                                                       15597, 15633, 37255,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66070, 0, 3,
                                                                       65035, 36625, 65080,
                                                                       15633, 15669, 37315,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66160, 0, 3,
                                                                       65080, 36655, 65125,
                                                                       15669, 15705, 37375,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66250, 0, 3,
                                                                       65125, 36685, 65170,
                                                                       15705, 15741, 37435,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66340, 0, 3,
                                                                       65170, 36715, 65215,
                                                                       15741, 15777, 37495,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66430, 0, 3,
                                                                       65215, 36745, 65260,
                                                                       15777, 15813, 37555,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66520, 0, 3,
                                                                       65260, 36775, 65305,
                                                                       15813, 15849, 37615,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 66610, 0, 3,
                                                                       65350, 36835, 65440,
                                                                       15921, 15981, 37675,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 66760, 0, 3,
                                                                       65440, 36895, 65530,
                                                                       15981, 16041, 37775,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 66910, 0, 3,
                                                                       65530, 36955, 65620,
                                                                       16041, 16101, 37875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67060, 0, 3,
                                                                       65620, 37015, 65710,
                                                                       16101, 16161, 37975,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67210, 0, 3,
                                                                       65710, 37075, 65800,
                                                                       16161, 16221, 38075,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67360, 0, 3,
                                                                       65800, 37135, 65890,
                                                                       16221, 16281, 38175,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67510, 0, 3,
                                                                       65890, 37195, 65980,
                                                                       16281, 16341, 38275,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67660, 0, 3,
                                                                       65980, 37255, 66070,
                                                                       16341, 16401, 38375,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67810, 0, 3,
                                                                       66070, 37315, 66160,
                                                                       16401, 16461, 38475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67960, 0, 3,
                                                                       66160, 37375, 66250,
                                                                       16461, 16521, 38575,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 68110, 0, 3,
                                                                       66250, 37435, 66340,
                                                                       16521, 16581, 38675,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 68260, 0, 3,
                                                                       66340, 37495, 66430,
                                                                       16581, 16641, 38775,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 68410, 0, 3,
                                                                       66430, 37555, 66520,
                                                                       16641, 16701, 38875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 68560, 0, 3,
                                                                       66610, 37675, 66760,
                                                                       16821, 16911, 38975,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 68785, 0, 3,
                                                                       66760, 37775, 66910,
                                                                       16911, 17001, 39125,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 69010, 0, 3,
                                                                       66910, 37875, 67060,
                                                                       17001, 17091, 39275,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 69235, 0, 3,
                                                                       67060, 37975, 67210,
                                                                       17091, 17181, 39425,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 69460, 0, 3,
                                                                       67210, 38075, 67360,
                                                                       17181, 17271, 39575,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 69685, 0, 3,
                                                                       67360, 38175, 67510,
                                                                       17271, 17361, 39725,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 69910, 0, 3,
                                                                       67510, 38275, 67660,
                                                                       17361, 17451, 39875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 70135, 0, 3,
                                                                       67660, 38375, 67810,
                                                                       17451, 17541, 40025,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 70360, 0, 3,
                                                                       67810, 38475, 67960,
                                                                       17541, 17631, 40175,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 70585, 0, 3,
                                                                       67960, 38575, 68110,
                                                                       17631, 17721, 40325,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 70810, 0, 3,
                                                                       68110, 38675, 68260,
                                                                       17721, 17811, 40475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 71035, 0, 3,
                                                                       68260, 38775, 68410,
                                                                       17811, 17901, 40625,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 71260, 0, 3,
                                                                       68560, 38975, 68785,
                                                                       18081, 18207, 40775,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 71575, 0, 3,
                                                                       68785, 39125, 69010,
                                                                       18207, 18333, 40985,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 71890, 0, 3,
                                                                       69010, 39275, 69235,
                                                                       18333, 18459, 41195,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 72205, 0, 3,
                                                                       69235, 39425, 69460,
                                                                       18459, 18585, 41405,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 72520, 0, 3,
                                                                       69460, 39575, 69685,
                                                                       18585, 18711, 41615,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 72835, 0, 3,
                                                                       69685, 39725, 69910,
                                                                       18711, 18837, 41825,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 73150, 0, 3,
                                                                       69910, 39875, 70135,
                                                                       18837, 18963, 42035,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 73465, 0, 3,
                                                                       70135, 40025, 70360,
                                                                       18963, 19089, 42245,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 73780, 0, 3,
                                                                       70360, 40175, 70585,
                                                                       19089, 19215, 42455,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 74095, 0, 3,
                                                                       70585, 40325, 70810,
                                                                       19215, 19341, 42665,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 74410, 0, 3,
                                                                       70810, 40475, 71035,
                                                                       19341, 19467, 42875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 74725, 0, 3,
                                                                       71260, 40775, 71575,
                                                                       19719, 19887, 43085,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 75145, 0, 3,
                                                                       71575, 40985, 71890,
                                                                       19887, 20055, 43365,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 75565, 0, 3,
                                                                       71890, 41195, 72205,
                                                                       20055, 20223, 43645,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 75985, 0, 3,
                                                                       72205, 41405, 72520,
                                                                       20223, 20391, 43925,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 76405, 0, 3,
                                                                       72520, 41615, 72835,
                                                                       20391, 20559, 44205,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 76825, 0, 3,
                                                                       72835, 41825, 73150,
                                                                       20559, 20727, 44485,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 77245, 0, 3,
                                                                       73150, 42035, 73465,
                                                                       20727, 20895, 44765,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 77665, 0, 3,
                                                                       73465, 42245, 73780,
                                                                       20895, 21063, 45045,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 78085, 0, 3,
                                                                       73780, 42455, 74095,
                                                                       21063, 21231, 45325,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 78505, 0, 3,
                                                                       74095, 42665, 74410,
                                                                       21231, 21399, 45605,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 78925, 0, 3,
                                                                       74725, 43085, 75145,
                                                                       21735, 21951, 45885,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 79465, 0, 3,
                                                                       75145, 43365, 75565,
                                                                       21951, 22167, 46245,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 80005, 0, 3,
                                                                       75565, 43645, 75985,
                                                                       22167, 22383, 46605,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 80545, 0, 3,
                                                                       75985, 43925, 76405,
                                                                       22383, 22599, 46965,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 81085, 0, 3,
                                                                       76405, 44205, 76825,
                                                                       22599, 22815, 47325,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 81625, 0, 3,
                                                                       76825, 44485, 77245,
                                                                       22815, 23031, 47685,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 82165, 0, 3,
                                                                       77245, 44765, 77665,
                                                                       23031, 23247, 48045,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 82705, 0, 3,
                                                                       77665, 45045, 78085,
                                                                       23247, 23463, 48405,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 83245, 0, 3,
                                                                       78085, 45325, 78505,
                                                                       23463, 23679, 48765,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 83785, 0, 3,
                                                                       78925, 45885, 79465,
                                                                       24111, 24381, 49125,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 84460, 0, 3,
                                                                       79465, 46245, 80005,
                                                                       24381, 24651, 49575,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 85135, 0, 3,
                                                                       80005, 46605, 80545,
                                                                       24651, 24921, 50025,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 85810, 0, 3,
                                                                       80545, 46965, 81085,
                                                                       24921, 25191, 50475,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 86485, 0, 3,
                                                                       81085, 47325, 81625,
                                                                       25191, 25461, 50925,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 87160, 0, 3,
                                                                       81625, 47685, 82165,
                                                                       25461, 25731, 51375,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 87835, 0, 3,
                                                                       82165, 48045, 82705,
                                                                       25731, 26001, 51825,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 88510, 0, 3,
                                                                       82705, 48405, 83245,
                                                                       26001, 26271, 52275,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 89185, 0, 3,
                                                                       83785, 49125, 84460,
                                                                       26811, 27141, 52725,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 90010, 0, 3,
                                                                       84460, 49575, 85135,
                                                                       27141, 27471, 53275,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 90835, 0, 3,
                                                                       85135, 50025, 85810,
                                                                       27471, 27801, 53825,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 91660, 0, 3,
                                                                       85810, 50475, 86485,
                                                                       27801, 28131, 54375,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 92485, 0, 3,
                                                                       86485, 50925, 87160,
                                                                       28131, 28461, 54925,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 93310, 0, 3,
                                                                       87160, 51375, 87835,
                                                                       28461, 28791, 55475,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 94135, 0, 3,
                                                                       87835, 51825, 88510,
                                                                       28791, 29121, 56025,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 94960, 0, 3,
                                                                       89185, 52725, 90010,
                                                                       29781, 30177, 56575,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 95950, 0, 3,
                                                                       90010, 53275, 90835,
                                                                       30177, 30573, 57235,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 96940, 0, 3,
                                                                       90835, 53825, 91660,
                                                                       30573, 30969, 57895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 97930, 0, 3,
                                                                       91660, 54375, 92485,
                                                                       30969, 31365, 58555,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 98920, 0, 3,
                                                                       92485, 54925, 93310,
                                                                       31365, 31761, 59215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 99910, 0, 3,
                                                                       93310, 55475, 94135,
                                                                       31761, 32157, 59875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 100900, 0, 3,
                                                                       94960, 56575, 95950,
                                                                       32949, 33417, 60535,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 102070, 0, 3,
                                                                       95950, 57235, 96940,
                                                                       33417, 33885, 61315,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 103240, 0, 3,
                                                                       96940, 57895, 97930,
                                                                       33885, 34353, 62095,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 104410, 0, 3,
                                                                       97930, 58555, 98920,
                                                                       34353, 34821, 62875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 105580, 0, 3,
                                                                       98920, 59215, 99910,
                                                                       34821, 35289, 63655,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106750, 3, 36225,
                                                                       36235, 64465, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106771, 3, 36235,
                                                                       36245, 64480, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106792, 3, 36245,
                                                                       36255, 64495, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106813, 3, 36255,
                                                                       36265, 64510, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106834, 3, 36265,
                                                                       36275, 64525, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106855, 3, 36275,
                                                                       36285, 64540, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106876, 3, 36285,
                                                                       36295, 64555, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106897, 3, 36295,
                                                                       36305, 64570, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106918, 3, 36305,
                                                                       36315, 64585, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106939, 3, 36315,
                                                                       36325, 64600, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106960, 3, 36325,
                                                                       36335, 64615, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106981, 3, 36335,
                                                                       36345, 64630, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 107002, 3, 36345,
                                                                       36355, 64645, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 107023, 3, 36355,
                                                                       36365, 64660, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107044, 0, 3,
                                                                       106750, 64465, 106771,
                                                                       36385, 36415, 64765,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107107, 0, 3,
                                                                       106771, 64480, 106792,
                                                                       36415, 36445, 64810,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107170, 0, 3,
                                                                       106792, 64495, 106813,
                                                                       36445, 36475, 64855,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107233, 0, 3,
                                                                       106813, 64510, 106834,
                                                                       36475, 36505, 64900,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107296, 0, 3,
                                                                       106834, 64525, 106855,
                                                                       36505, 36535, 64945,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107359, 0, 3,
                                                                       106855, 64540, 106876,
                                                                       36535, 36565, 64990,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107422, 0, 3,
                                                                       106876, 64555, 106897,
                                                                       36565, 36595, 65035,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107485, 0, 3,
                                                                       106897, 64570, 106918,
                                                                       36595, 36625, 65080,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107548, 0, 3,
                                                                       106918, 64585, 106939,
                                                                       36625, 36655, 65125,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107611, 0, 3,
                                                                       106939, 64600, 106960,
                                                                       36655, 36685, 65170,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107674, 0, 3,
                                                                       106960, 64615, 106981,
                                                                       36685, 36715, 65215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107737, 0, 3,
                                                                       106981, 64630, 107002,
                                                                       36715, 36745, 65260,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107800, 0, 3,
                                                                       107002, 64645, 107023,
                                                                       36745, 36775, 65305,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 107863, 0, 3,
                                                                       107044, 64765, 107107,
                                                                       36835, 36895, 65530,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 107989, 0, 3,
                                                                       107107, 64810, 107170,
                                                                       36895, 36955, 65620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108115, 0, 3,
                                                                       107170, 64855, 107233,
                                                                       36955, 37015, 65710,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108241, 0, 3,
                                                                       107233, 64900, 107296,
                                                                       37015, 37075, 65800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108367, 0, 3,
                                                                       107296, 64945, 107359,
                                                                       37075, 37135, 65890,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108493, 0, 3,
                                                                       107359, 64990, 107422,
                                                                       37135, 37195, 65980,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108619, 0, 3,
                                                                       107422, 65035, 107485,
                                                                       37195, 37255, 66070,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108745, 0, 3,
                                                                       107485, 65080, 107548,
                                                                       37255, 37315, 66160,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108871, 0, 3,
                                                                       107548, 65125, 107611,
                                                                       37315, 37375, 66250,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108997, 0, 3,
                                                                       107611, 65170, 107674,
                                                                       37375, 37435, 66340,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 109123, 0, 3,
                                                                       107674, 65215, 107737,
                                                                       37435, 37495, 66430,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 109249, 0, 3,
                                                                       107737, 65260, 107800,
                                                                       37495, 37555, 66520,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109375, 0, 3,
                                                                       107863, 65530, 107989,
                                                                       37675, 37775, 66910,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109585, 0, 3,
                                                                       107989, 65620, 108115,
                                                                       37775, 37875, 67060,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109795, 0, 3,
                                                                       108115, 65710, 108241,
                                                                       37875, 37975, 67210,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 110005, 0, 3,
                                                                       108241, 65800, 108367,
                                                                       37975, 38075, 67360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 110215, 0, 3,
                                                                       108367, 65890, 108493,
                                                                       38075, 38175, 67510,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 110425, 0, 3,
                                                                       108493, 65980, 108619,
                                                                       38175, 38275, 67660,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 110635, 0, 3,
                                                                       108619, 66070, 108745,
                                                                       38275, 38375, 67810,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 110845, 0, 3,
                                                                       108745, 66160, 108871,
                                                                       38375, 38475, 67960,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 111055, 0, 3,
                                                                       108871, 66250, 108997,
                                                                       38475, 38575, 68110,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 111265, 0, 3,
                                                                       108997, 66340, 109123,
                                                                       38575, 38675, 68260,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 111475, 0, 3,
                                                                       109123, 66430, 109249,
                                                                       38675, 38775, 68410,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 111685, 0, 3,
                                                                       109375, 66910, 109585,
                                                                       38975, 39125, 69010,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112000, 0, 3,
                                                                       109585, 67060, 109795,
                                                                       39125, 39275, 69235,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112315, 0, 3,
                                                                       109795, 67210, 110005,
                                                                       39275, 39425, 69460,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112630, 0, 3,
                                                                       110005, 67360, 110215,
                                                                       39425, 39575, 69685,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112945, 0, 3,
                                                                       110215, 67510, 110425,
                                                                       39575, 39725, 69910,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 113260, 0, 3,
                                                                       110425, 67660, 110635,
                                                                       39725, 39875, 70135,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 113575, 0, 3,
                                                                       110635, 67810, 110845,
                                                                       39875, 40025, 70360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 113890, 0, 3,
                                                                       110845, 67960, 111055,
                                                                       40025, 40175, 70585,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 114205, 0, 3,
                                                                       111055, 68110, 111265,
                                                                       40175, 40325, 70810,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 114520, 0, 3,
                                                                       111265, 68260, 111475,
                                                                       40325, 40475, 71035,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 114835, 0, 3,
                                                                       111685, 69010, 112000,
                                                                       40775, 40985, 71890,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 115276, 0, 3,
                                                                       112000, 69235, 112315,
                                                                       40985, 41195, 72205,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 115717, 0, 3,
                                                                       112315, 69460, 112630,
                                                                       41195, 41405, 72520,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 116158, 0, 3,
                                                                       112630, 69685, 112945,
                                                                       41405, 41615, 72835,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 116599, 0, 3,
                                                                       112945, 69910, 113260,
                                                                       41615, 41825, 73150,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 117040, 0, 3,
                                                                       113260, 70135, 113575,
                                                                       41825, 42035, 73465,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 117481, 0, 3,
                                                                       113575, 70360, 113890,
                                                                       42035, 42245, 73780,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 117922, 0, 3,
                                                                       113890, 70585, 114205,
                                                                       42245, 42455, 74095,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 118363, 0, 3,
                                                                       114205, 70810, 114520,
                                                                       42455, 42665, 74410,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 118804, 0, 3,
                                                                       114835, 71890, 115276,
                                                                       43085, 43365, 75565,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 119392, 0, 3,
                                                                       115276, 72205, 115717,
                                                                       43365, 43645, 75985,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 119980, 0, 3,
                                                                       115717, 72520, 116158,
                                                                       43645, 43925, 76405,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 120568, 0, 3,
                                                                       116158, 72835, 116599,
                                                                       43925, 44205, 76825,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 121156, 0, 3,
                                                                       116599, 73150, 117040,
                                                                       44205, 44485, 77245,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 121744, 0, 3,
                                                                       117040, 73465, 117481,
                                                                       44485, 44765, 77665,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 122332, 0, 3,
                                                                       117481, 73780, 117922,
                                                                       44765, 45045, 78085,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 122920, 0, 3,
                                                                       117922, 74095, 118363,
                                                                       45045, 45325, 78505,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 123508, 0, 3,
                                                                       118804, 75565, 119392,
                                                                       45885, 46245, 80005,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 124264, 0, 3,
                                                                       119392, 75985, 119980,
                                                                       46245, 46605, 80545,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 125020, 0, 3,
                                                                       119980, 76405, 120568,
                                                                       46605, 46965, 81085,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 125776, 0, 3,
                                                                       120568, 76825, 121156,
                                                                       46965, 47325, 81625,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 126532, 0, 3,
                                                                       121156, 77245, 121744,
                                                                       47325, 47685, 82165,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 127288, 0, 3,
                                                                       121744, 77665, 122332,
                                                                       47685, 48045, 82705,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 128044, 0, 3,
                                                                       122332, 78085, 122920,
                                                                       48045, 48405, 83245,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 128800, 0, 3,
                                                                       123508, 80005, 124264,
                                                                       49125, 49575, 85135,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 129745, 0, 3,
                                                                       124264, 80545, 125020,
                                                                       49575, 50025, 85810,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 130690, 0, 3,
                                                                       125020, 81085, 125776,
                                                                       50025, 50475, 86485,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 131635, 0, 3,
                                                                       125776, 81625, 126532,
                                                                       50475, 50925, 87160,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 132580, 0, 3,
                                                                       126532, 82165, 127288,
                                                                       50925, 51375, 87835,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 133525, 0, 3,
                                                                       127288, 82705, 128044,
                                                                       51375, 51825, 88510,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 134470, 0, 3,
                                                                       128800, 85135, 129745,
                                                                       52725, 53275, 90835,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 135625, 0, 3,
                                                                       129745, 85810, 130690,
                                                                       53275, 53825, 91660,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 136780, 0, 3,
                                                                       130690, 86485, 131635,
                                                                       53825, 54375, 92485,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 137935, 0, 3,
                                                                       131635, 87160, 132580,
                                                                       54375, 54925, 93310,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 139090, 0, 3,
                                                                       132580, 87835, 133525,
                                                                       54925, 55475, 94135,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 140245, 0, 3,
                                                                       134470, 90835, 135625,
                                                                       56575, 57235, 96940,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 141631, 0, 3,
                                                                       135625, 91660, 136780,
                                                                       57235, 57895, 97930,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 143017, 0, 3,
                                                                       136780, 92485, 137935,
                                                                       57895, 58555, 98920,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 144403, 0, 3,
                                                                       137935, 93310, 139090,
                                                                       58555, 59215, 99910,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 145789, 0, 3,
                                                                       140245, 96940, 141631,
                                                                       60535, 61315, 103240,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 147427, 0, 3,
                                                                       141631, 97930, 143017,
                                                                       61315, 62095, 104410,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 149065, 0, 3,
                                                                       143017, 98920, 144403,
                                                                       62095, 62875, 105580,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150703, 3, 64435,
                                                                       64450, 106750, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150731, 3, 64450,
                                                                       64465, 106771, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150759, 3, 64465,
                                                                       64480, 106792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150787, 3, 64480,
                                                                       64495, 106813, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150815, 3, 64495,
                                                                       64510, 106834, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150843, 3, 64510,
                                                                       64525, 106855, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150871, 3, 64525,
                                                                       64540, 106876, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150899, 3, 64540,
                                                                       64555, 106897, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150927, 3, 64555,
                                                                       64570, 106918, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150955, 3, 64570,
                                                                       64585, 106939, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150983, 3, 64585,
                                                                       64600, 106960, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 151011, 3, 64600,
                                                                       64615, 106981, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 151039, 3, 64615,
                                                                       64630, 107002, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 151067, 3, 64630,
                                                                       64645, 107023, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151095, 0, 3,
                                                                       150703, 106750, 150731,
                                                                       64675, 64720, 107044,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151179, 0, 3,
                                                                       150731, 106771, 150759,
                                                                       64720, 64765, 107107,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151263, 0, 3,
                                                                       150759, 106792, 150787,
                                                                       64765, 64810, 107170,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151347, 0, 3,
                                                                       150787, 106813, 150815,
                                                                       64810, 64855, 107233,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151431, 0, 3,
                                                                       150815, 106834, 150843,
                                                                       64855, 64900, 107296,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151515, 0, 3,
                                                                       150843, 106855, 150871,
                                                                       64900, 64945, 107359,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151599, 0, 3,
                                                                       150871, 106876, 150899,
                                                                       64945, 64990, 107422,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151683, 0, 3,
                                                                       150899, 106897, 150927,
                                                                       64990, 65035, 107485,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151767, 0, 3,
                                                                       150927, 106918, 150955,
                                                                       65035, 65080, 107548,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151851, 0, 3,
                                                                       150955, 106939, 150983,
                                                                       65080, 65125, 107611,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151935, 0, 3,
                                                                       150983, 106960, 151011,
                                                                       65125, 65170, 107674,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 152019, 0, 3,
                                                                       151011, 106981, 151039,
                                                                       65170, 65215, 107737,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 152103, 0, 3,
                                                                       151039, 107002, 151067,
                                                                       65215, 65260, 107800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152187, 0, 3,
                                                                       151095, 107044, 151179,
                                                                       65350, 65440, 107863,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152355, 0, 3,
                                                                       151179, 107107, 151263,
                                                                       65440, 65530, 107989,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152523, 0, 3,
                                                                       151263, 107170, 151347,
                                                                       65530, 65620, 108115,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152691, 0, 3,
                                                                       151347, 107233, 151431,
                                                                       65620, 65710, 108241,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152859, 0, 3,
                                                                       151431, 107296, 151515,
                                                                       65710, 65800, 108367,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153027, 0, 3,
                                                                       151515, 107359, 151599,
                                                                       65800, 65890, 108493,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153195, 0, 3,
                                                                       151599, 107422, 151683,
                                                                       65890, 65980, 108619,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153363, 0, 3,
                                                                       151683, 107485, 151767,
                                                                       65980, 66070, 108745,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153531, 0, 3,
                                                                       151767, 107548, 151851,
                                                                       66070, 66160, 108871,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153699, 0, 3,
                                                                       151851, 107611, 151935,
                                                                       66160, 66250, 108997,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153867, 0, 3,
                                                                       151935, 107674, 152019,
                                                                       66250, 66340, 109123,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 154035, 0, 3,
                                                                       152019, 107737, 152103,
                                                                       66340, 66430, 109249,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 154203, 0, 3,
                                                                       152187, 107863, 152355,
                                                                       66610, 66760, 109375,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 154483, 0, 3,
                                                                       152355, 107989, 152523,
                                                                       66760, 66910, 109585,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 154763, 0, 3,
                                                                       152523, 108115, 152691,
                                                                       66910, 67060, 109795,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155043, 0, 3,
                                                                       152691, 108241, 152859,
                                                                       67060, 67210, 110005,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155323, 0, 3,
                                                                       152859, 108367, 153027,
                                                                       67210, 67360, 110215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155603, 0, 3,
                                                                       153027, 108493, 153195,
                                                                       67360, 67510, 110425,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155883, 0, 3,
                                                                       153195, 108619, 153363,
                                                                       67510, 67660, 110635,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 156163, 0, 3,
                                                                       153363, 108745, 153531,
                                                                       67660, 67810, 110845,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 156443, 0, 3,
                                                                       153531, 108871, 153699,
                                                                       67810, 67960, 111055,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 156723, 0, 3,
                                                                       153699, 108997, 153867,
                                                                       67960, 68110, 111265,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 157003, 0, 3,
                                                                       153867, 109123, 154035,
                                                                       68110, 68260, 111475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 157283, 0, 3,
                                                                       154203, 109375, 154483,
                                                                       68560, 68785, 111685,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 157703, 0, 3,
                                                                       154483, 109585, 154763,
                                                                       68785, 69010, 112000,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 158123, 0, 3,
                                                                       154763, 109795, 155043,
                                                                       69010, 69235, 112315,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 158543, 0, 3,
                                                                       155043, 110005, 155323,
                                                                       69235, 69460, 112630,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 158963, 0, 3,
                                                                       155323, 110215, 155603,
                                                                       69460, 69685, 112945,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 159383, 0, 3,
                                                                       155603, 110425, 155883,
                                                                       69685, 69910, 113260,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 159803, 0, 3,
                                                                       155883, 110635, 156163,
                                                                       69910, 70135, 113575,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 160223, 0, 3,
                                                                       156163, 110845, 156443,
                                                                       70135, 70360, 113890,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 160643, 0, 3,
                                                                       156443, 111055, 156723,
                                                                       70360, 70585, 114205,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 161063, 0, 3,
                                                                       156723, 111265, 157003,
                                                                       70585, 70810, 114520,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 161483, 0, 3,
                                                                       157283, 111685, 157703,
                                                                       71260, 71575, 114835,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 162071, 0, 3,
                                                                       157703, 112000, 158123,
                                                                       71575, 71890, 115276,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 162659, 0, 3,
                                                                       158123, 112315, 158543,
                                                                       71890, 72205, 115717,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 163247, 0, 3,
                                                                       158543, 112630, 158963,
                                                                       72205, 72520, 116158,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 163835, 0, 3,
                                                                       158963, 112945, 159383,
                                                                       72520, 72835, 116599,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 164423, 0, 3,
                                                                       159383, 113260, 159803,
                                                                       72835, 73150, 117040,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 165011, 0, 3,
                                                                       159803, 113575, 160223,
                                                                       73150, 73465, 117481,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 165599, 0, 3,
                                                                       160223, 113890, 160643,
                                                                       73465, 73780, 117922,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 166187, 0, 3,
                                                                       160643, 114205, 161063,
                                                                       73780, 74095, 118363,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 166775, 0, 3,
                                                                       161483, 114835, 162071,
                                                                       74725, 75145, 118804,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 167559, 0, 3,
                                                                       162071, 115276, 162659,
                                                                       75145, 75565, 119392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 168343, 0, 3,
                                                                       162659, 115717, 163247,
                                                                       75565, 75985, 119980,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 169127, 0, 3,
                                                                       163247, 116158, 163835,
                                                                       75985, 76405, 120568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 169911, 0, 3,
                                                                       163835, 116599, 164423,
                                                                       76405, 76825, 121156,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 170695, 0, 3,
                                                                       164423, 117040, 165011,
                                                                       76825, 77245, 121744,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 171479, 0, 3,
                                                                       165011, 117481, 165599,
                                                                       77245, 77665, 122332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 172263, 0, 3,
                                                                       165599, 117922, 166187,
                                                                       77665, 78085, 122920,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 173047, 0, 3,
                                                                       166775, 118804, 167559,
                                                                       78925, 79465, 123508,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 174055, 0, 3,
                                                                       167559, 119392, 168343,
                                                                       79465, 80005, 124264,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 175063, 0, 3,
                                                                       168343, 119980, 169127,
                                                                       80005, 80545, 125020,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 176071, 0, 3,
                                                                       169127, 120568, 169911,
                                                                       80545, 81085, 125776,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 177079, 0, 3,
                                                                       169911, 121156, 170695,
                                                                       81085, 81625, 126532,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 178087, 0, 3,
                                                                       170695, 121744, 171479,
                                                                       81625, 82165, 127288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 179095, 0, 3,
                                                                       171479, 122332, 172263,
                                                                       82165, 82705, 128044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 180103, 0, 3,
                                                                       173047, 123508, 174055,
                                                                       83785, 84460, 128800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 181363, 0, 3,
                                                                       174055, 124264, 175063,
                                                                       84460, 85135, 129745,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 182623, 0, 3,
                                                                       175063, 125020, 176071,
                                                                       85135, 85810, 130690,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 183883, 0, 3,
                                                                       176071, 125776, 177079,
                                                                       85810, 86485, 131635,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 185143, 0, 3,
                                                                       177079, 126532, 178087,
                                                                       86485, 87160, 132580,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 186403, 0, 3,
                                                                       178087, 127288, 179095,
                                                                       87160, 87835, 133525,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 187663, 0, 3,
                                                                       180103, 128800, 181363,
                                                                       89185, 90010, 134470,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 189203, 0, 3,
                                                                       181363, 129745, 182623,
                                                                       90010, 90835, 135625,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 190743, 0, 3,
                                                                       182623, 130690, 183883,
                                                                       90835, 91660, 136780,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 192283, 0, 3,
                                                                       183883, 131635, 185143,
                                                                       91660, 92485, 137935,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 193823, 0, 3,
                                                                       185143, 132580, 186403,
                                                                       92485, 93310, 139090,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 195363, 0, 3,
                                                                       187663, 134470, 189203,
                                                                       94960, 95950, 140245,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 197211, 0, 3,
                                                                       189203, 135625, 190743,
                                                                       95950, 96940, 141631,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 199059, 0, 3,
                                                                       190743, 136780, 192283,
                                                                       96940, 97930, 143017,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 200907, 0, 3,
                                                                       192283, 137935, 193823,
                                                                       97930, 98920, 144403,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 202755, 0, 3,
                                                                       195363, 140245, 197211,
                                                                       100900, 102070, 145789,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 204939, 0, 3,
                                                                       197211, 141631, 199059,
                                                                       102070, 103240, 147427,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 207123, 0, 3,
                                                                       199059, 143017, 200907,
                                                                       103240, 104410, 149065,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209307, 3, 106750,
                                                                       106771, 150759, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209343, 3, 106771,
                                                                       106792, 150787, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209379, 3, 106792,
                                                                       106813, 150815, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209415, 3, 106813,
                                                                       106834, 150843, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209451, 3, 106834,
                                                                       106855, 150871, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209487, 3, 106855,
                                                                       106876, 150899, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209523, 3, 106876,
                                                                       106897, 150927, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209559, 3, 106897,
                                                                       106918, 150955, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209595, 3, 106918,
                                                                       106939, 150983, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209631, 3, 106939,
                                                                       106960, 151011, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209667, 3, 106960,
                                                                       106981, 151039, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209703, 3, 106981,
                                                                       107002, 151067, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 209739, 0, 3,
                                                                       209307, 150759, 209343,
                                                                       107044, 107107, 151263,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 209847, 0, 3,
                                                                       209343, 150787, 209379,
                                                                       107107, 107170, 151347,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 209955, 0, 3,
                                                                       209379, 150815, 209415,
                                                                       107170, 107233, 151431,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210063, 0, 3,
                                                                       209415, 150843, 209451,
                                                                       107233, 107296, 151515,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210171, 0, 3,
                                                                       209451, 150871, 209487,
                                                                       107296, 107359, 151599,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210279, 0, 3,
                                                                       209487, 150899, 209523,
                                                                       107359, 107422, 151683,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210387, 0, 3,
                                                                       209523, 150927, 209559,
                                                                       107422, 107485, 151767,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210495, 0, 3,
                                                                       209559, 150955, 209595,
                                                                       107485, 107548, 151851,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210603, 0, 3,
                                                                       209595, 150983, 209631,
                                                                       107548, 107611, 151935,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210711, 0, 3,
                                                                       209631, 151011, 209667,
                                                                       107611, 107674, 152019,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210819, 0, 3,
                                                                       209667, 151039, 209703,
                                                                       107674, 107737, 152103,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 210927, 0, 3,
                                                                       209739, 151263, 209847,
                                                                       107863, 107989, 152523,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 211143, 0, 3,
                                                                       209847, 151347, 209955,
                                                                       107989, 108115, 152691,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 211359, 0, 3,
                                                                       209955, 151431, 210063,
                                                                       108115, 108241, 152859,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 211575, 0, 3,
                                                                       210063, 151515, 210171,
                                                                       108241, 108367, 153027,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 211791, 0, 3,
                                                                       210171, 151599, 210279,
                                                                       108367, 108493, 153195,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 212007, 0, 3,
                                                                       210279, 151683, 210387,
                                                                       108493, 108619, 153363,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 212223, 0, 3,
                                                                       210387, 151767, 210495,
                                                                       108619, 108745, 153531,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 212439, 0, 3,
                                                                       210495, 151851, 210603,
                                                                       108745, 108871, 153699,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 212655, 0, 3,
                                                                       210603, 151935, 210711,
                                                                       108871, 108997, 153867,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 212871, 0, 3,
                                                                       210711, 152019, 210819,
                                                                       108997, 109123, 154035,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 213087, 0, 3,
                                                                       210927, 152523, 211143,
                                                                       109375, 109585, 154763,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 213447, 0, 3,
                                                                       211143, 152691, 211359,
                                                                       109585, 109795, 155043,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 213807, 0, 3,
                                                                       211359, 152859, 211575,
                                                                       109795, 110005, 155323,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 214167, 0, 3,
                                                                       211575, 153027, 211791,
                                                                       110005, 110215, 155603,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 214527, 0, 3,
                                                                       211791, 153195, 212007,
                                                                       110215, 110425, 155883,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 214887, 0, 3,
                                                                       212007, 153363, 212223,
                                                                       110425, 110635, 156163,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 215247, 0, 3,
                                                                       212223, 153531, 212439,
                                                                       110635, 110845, 156443,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 215607, 0, 3,
                                                                       212439, 153699, 212655,
                                                                       110845, 111055, 156723,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 215967, 0, 3,
                                                                       212655, 153867, 212871,
                                                                       111055, 111265, 157003,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 216327, 0, 3,
                                                                       213087, 154763, 213447,
                                                                       111685, 112000, 158123,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 216867, 0, 3,
                                                                       213447, 155043, 213807,
                                                                       112000, 112315, 158543,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 217407, 0, 3,
                                                                       213807, 155323, 214167,
                                                                       112315, 112630, 158963,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 217947, 0, 3,
                                                                       214167, 155603, 214527,
                                                                       112630, 112945, 159383,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 218487, 0, 3,
                                                                       214527, 155883, 214887,
                                                                       112945, 113260, 159803,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 219027, 0, 3,
                                                                       214887, 156163, 215247,
                                                                       113260, 113575, 160223,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 219567, 0, 3,
                                                                       215247, 156443, 215607,
                                                                       113575, 113890, 160643,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 220107, 0, 3,
                                                                       215607, 156723, 215967,
                                                                       113890, 114205, 161063,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 220647, 0, 3,
                                                                       216327, 158123, 216867,
                                                                       114835, 115276, 162659,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 221403, 0, 3,
                                                                       216867, 158543, 217407,
                                                                       115276, 115717, 163247,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 222159, 0, 3,
                                                                       217407, 158963, 217947,
                                                                       115717, 116158, 163835,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 222915, 0, 3,
                                                                       217947, 159383, 218487,
                                                                       116158, 116599, 164423,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 223671, 0, 3,
                                                                       218487, 159803, 219027,
                                                                       116599, 117040, 165011,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 224427, 0, 3,
                                                                       219027, 160223, 219567,
                                                                       117040, 117481, 165599,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 225183, 0, 3,
                                                                       219567, 160643, 220107,
                                                                       117481, 117922, 166187,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 225939, 0, 3,
                                                                       220647, 162659, 221403,
                                                                       118804, 119392, 168343,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 226947, 0, 3,
                                                                       221403, 163247, 222159,
                                                                       119392, 119980, 169127,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 227955, 0, 3,
                                                                       222159, 163835, 222915,
                                                                       119980, 120568, 169911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 228963, 0, 3,
                                                                       222915, 164423, 223671,
                                                                       120568, 121156, 170695,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 229971, 0, 3,
                                                                       223671, 165011, 224427,
                                                                       121156, 121744, 171479,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 230979, 0, 3,
                                                                       224427, 165599, 225183,
                                                                       121744, 122332, 172263,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 231987, 0, 3,
                                                                       225939, 168343, 226947,
                                                                       123508, 124264, 175063,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 233283, 0, 3,
                                                                       226947, 169127, 227955,
                                                                       124264, 125020, 176071,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 234579, 0, 3,
                                                                       227955, 169911, 228963,
                                                                       125020, 125776, 177079,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 235875, 0, 3,
                                                                       228963, 170695, 229971,
                                                                       125776, 126532, 178087,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 237171, 0, 3,
                                                                       229971, 171479, 230979,
                                                                       126532, 127288, 179095,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 238467, 0, 3,
                                                                       231987, 175063, 233283,
                                                                       128800, 129745, 182623,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 240087, 0, 3,
                                                                       233283, 176071, 234579,
                                                                       129745, 130690, 183883,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 241707, 0, 3,
                                                                       234579, 177079, 235875,
                                                                       130690, 131635, 185143,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 243327, 0, 3,
                                                                       235875, 178087, 237171,
                                                                       131635, 132580, 186403,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 244947, 0, 3,
                                                                       238467, 182623, 240087,
                                                                       134470, 135625, 190743,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 246927, 0, 3,
                                                                       240087, 183883, 241707,
                                                                       135625, 136780, 192283,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 248907, 0, 3,
                                                                       241707, 185143, 243327,
                                                                       136780, 137935, 193823,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 250887, 0, 3,
                                                                       244947, 190743, 246927,
                                                                       140245, 141631, 199059,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 253263, 0, 3,
                                                                       246927, 192283, 248907,
                                                                       141631, 143017, 200907,
                                                                       ncols, gamma, p, q);

                    compute_prim_sok_three_center_electron_repulsion_0(buffer, 255639, 0, 3,
                                                                       250887, 199059, 253263,
                                                                       145789, 147427, 207123,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258447, 3, 150703,
                                                                       150731, 209307, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258492, 3, 150731,
                                                                       150759, 209343, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258537, 3, 150759,
                                                                       150787, 209379, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258582, 3, 150787,
                                                                       150815, 209415, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258627, 3, 150815,
                                                                       150843, 209451, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258672, 3, 150843,
                                                                       150871, 209487, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258717, 3, 150871,
                                                                       150899, 209523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258762, 3, 150899,
                                                                       150927, 209559, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258807, 3, 150927,
                                                                       150955, 209595, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258852, 3, 150955,
                                                                       150983, 209631, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258897, 3, 150983,
                                                                       151011, 209667, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258942, 3, 151011,
                                                                       151039, 209703, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 258987, 0, 3,
                                                                       258447, 209307, 258492,
                                                                       151095, 151179, 209739,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259122, 0, 3,
                                                                       258492, 209343, 258537,
                                                                       151179, 151263, 209847,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259257, 0, 3,
                                                                       258537, 209379, 258582,
                                                                       151263, 151347, 209955,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259392, 0, 3,
                                                                       258582, 209415, 258627,
                                                                       151347, 151431, 210063,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259527, 0, 3,
                                                                       258627, 209451, 258672,
                                                                       151431, 151515, 210171,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259662, 0, 3,
                                                                       258672, 209487, 258717,
                                                                       151515, 151599, 210279,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259797, 0, 3,
                                                                       258717, 209523, 258762,
                                                                       151599, 151683, 210387,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259932, 0, 3,
                                                                       258762, 209559, 258807,
                                                                       151683, 151767, 210495,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 260067, 0, 3,
                                                                       258807, 209595, 258852,
                                                                       151767, 151851, 210603,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 260202, 0, 3,
                                                                       258852, 209631, 258897,
                                                                       151851, 151935, 210711,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 260337, 0, 3,
                                                                       258897, 209667, 258942,
                                                                       151935, 152019, 210819,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 260472, 0, 3,
                                                                       258987, 209739, 259122,
                                                                       152187, 152355, 210927,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 260742, 0, 3,
                                                                       259122, 209847, 259257,
                                                                       152355, 152523, 211143,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 261012, 0, 3,
                                                                       259257, 209955, 259392,
                                                                       152523, 152691, 211359,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 261282, 0, 3,
                                                                       259392, 210063, 259527,
                                                                       152691, 152859, 211575,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 261552, 0, 3,
                                                                       259527, 210171, 259662,
                                                                       152859, 153027, 211791,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 261822, 0, 3,
                                                                       259662, 210279, 259797,
                                                                       153027, 153195, 212007,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 262092, 0, 3,
                                                                       259797, 210387, 259932,
                                                                       153195, 153363, 212223,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 262362, 0, 3,
                                                                       259932, 210495, 260067,
                                                                       153363, 153531, 212439,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 262632, 0, 3,
                                                                       260067, 210603, 260202,
                                                                       153531, 153699, 212655,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 262902, 0, 3,
                                                                       260202, 210711, 260337,
                                                                       153699, 153867, 212871,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 263172, 0, 3,
                                                                       260472, 210927, 260742,
                                                                       154203, 154483, 213087,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 263622, 0, 3,
                                                                       260742, 211143, 261012,
                                                                       154483, 154763, 213447,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 264072, 0, 3,
                                                                       261012, 211359, 261282,
                                                                       154763, 155043, 213807,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 264522, 0, 3,
                                                                       261282, 211575, 261552,
                                                                       155043, 155323, 214167,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 264972, 0, 3,
                                                                       261552, 211791, 261822,
                                                                       155323, 155603, 214527,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 265422, 0, 3,
                                                                       261822, 212007, 262092,
                                                                       155603, 155883, 214887,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 265872, 0, 3,
                                                                       262092, 212223, 262362,
                                                                       155883, 156163, 215247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 266322, 0, 3,
                                                                       262362, 212439, 262632,
                                                                       156163, 156443, 215607,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 266772, 0, 3,
                                                                       262632, 212655, 262902,
                                                                       156443, 156723, 215967,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 267222, 0, 3,
                                                                       263172, 213087, 263622,
                                                                       157283, 157703, 216327,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 267897, 0, 3,
                                                                       263622, 213447, 264072,
                                                                       157703, 158123, 216867,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 268572, 0, 3,
                                                                       264072, 213807, 264522,
                                                                       158123, 158543, 217407,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 269247, 0, 3,
                                                                       264522, 214167, 264972,
                                                                       158543, 158963, 217947,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 269922, 0, 3,
                                                                       264972, 214527, 265422,
                                                                       158963, 159383, 218487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 270597, 0, 3,
                                                                       265422, 214887, 265872,
                                                                       159383, 159803, 219027,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 271272, 0, 3,
                                                                       265872, 215247, 266322,
                                                                       159803, 160223, 219567,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 271947, 0, 3,
                                                                       266322, 215607, 266772,
                                                                       160223, 160643, 220107,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 272622, 0, 3,
                                                                       267222, 216327, 267897,
                                                                       161483, 162071, 220647,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 273567, 0, 3,
                                                                       267897, 216867, 268572,
                                                                       162071, 162659, 221403,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 274512, 0, 3,
                                                                       268572, 217407, 269247,
                                                                       162659, 163247, 222159,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 275457, 0, 3,
                                                                       269247, 217947, 269922,
                                                                       163247, 163835, 222915,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 276402, 0, 3,
                                                                       269922, 218487, 270597,
                                                                       163835, 164423, 223671,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 277347, 0, 3,
                                                                       270597, 219027, 271272,
                                                                       164423, 165011, 224427,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 278292, 0, 3,
                                                                       271272, 219567, 271947,
                                                                       165011, 165599, 225183,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 279237, 0, 3,
                                                                       272622, 220647, 273567,
                                                                       166775, 167559, 225939,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 280497, 0, 3,
                                                                       273567, 221403, 274512,
                                                                       167559, 168343, 226947,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 281757, 0, 3,
                                                                       274512, 222159, 275457,
                                                                       168343, 169127, 227955,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 283017, 0, 3,
                                                                       275457, 222915, 276402,
                                                                       169127, 169911, 228963,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 284277, 0, 3,
                                                                       276402, 223671, 277347,
                                                                       169911, 170695, 229971,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 285537, 0, 3,
                                                                       277347, 224427, 278292,
                                                                       170695, 171479, 230979,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 286797, 0, 3,
                                                                       279237, 225939, 280497,
                                                                       173047, 174055, 231987,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 288417, 0, 3,
                                                                       280497, 226947, 281757,
                                                                       174055, 175063, 233283,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 290037, 0, 3,
                                                                       281757, 227955, 283017,
                                                                       175063, 176071, 234579,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 291657, 0, 3,
                                                                       283017, 228963, 284277,
                                                                       176071, 177079, 235875,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 293277, 0, 3,
                                                                       284277, 229971, 285537,
                                                                       177079, 178087, 237171,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 294897, 0, 3,
                                                                       286797, 231987, 288417,
                                                                       180103, 181363, 238467,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 296922, 0, 3,
                                                                       288417, 233283, 290037,
                                                                       181363, 182623, 240087,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 298947, 0, 3,
                                                                       290037, 234579, 291657,
                                                                       182623, 183883, 241707,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 300972, 0, 3,
                                                                       291657, 235875, 293277,
                                                                       183883, 185143, 243327,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 302997, 0, 3,
                                                                       294897, 238467, 296922,
                                                                       187663, 189203, 244947,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 305472, 0, 3,
                                                                       296922, 240087, 298947,
                                                                       189203, 190743, 246927,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 307947, 0, 3,
                                                                       298947, 241707, 300972,
                                                                       190743, 192283, 248907,
                                                                       ncols, gamma, p, q);

                    compute_prim_snl_three_center_electron_repulsion_0(buffer, 310422, 0, 3,
                                                                       302997, 244947, 305472,
                                                                       195363, 197211, 250887,
                                                                       ncols, gamma, p, q);

                    compute_prim_snl_three_center_electron_repulsion_0(buffer, 313392, 0, 3,
                                                                       305472, 246927, 307947,
                                                                       197211, 199059, 253263,
                                                                       ncols, gamma, p, q);

                    compute_prim_sol_three_center_electron_repulsion_0(buffer, 316362, 0, 3,
                                                                       310422, 250887, 313392,
                                                                       202755, 204939, 255639,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 319872, 279237, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 321608, 286797, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 323840, 294897, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 326630, 302997, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 330040, 310422, 2970, ncols);

                    simdfunc::contract_primitives(buffer, 334132, 316362, 3510, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 321132, 319872, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 323228, 321608, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 325865, 323840, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 329105, 326630, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 333010, 330040, 66, 1, nmax);

        simdtrf::transform_l_inner(buffer, 337642, 334132, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 338968, 321132, 323228, 17, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 340396, 323228, 325865, 17, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 342232, 325865, 329105, 17, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 344527, 329105, 333010, 17, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 347332, 333010, 337642, 17, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 350698, 338968, 340396, 17, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 353554, 340396, 342232, 17, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 357226, 342232, 344527, 17, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 361816, 344527, 347332, 17, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 367426, 350698, 353554, 17, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 372186, 353554, 357226, 17, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 378306, 357226, 361816, 17, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 385956, 367426, 372186, 17, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 393096, 372186, 378306, 17, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 402276, 385956, 393096, 17, nmax);

        simdtrf::transform_i_inner(buffer, 412272, 402276, 21, 17, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 412272, 221, nmax);
    }

    for (size_t m = 0; m < 2431; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
