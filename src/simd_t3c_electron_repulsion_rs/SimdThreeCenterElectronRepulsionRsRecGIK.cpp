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


#include "SimdThreeCenterElectronRepulsionRsRecGIK.hpp"

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
compute_rs_gik_three_center_electron_repulsion(double               *values,
                                               const size_t          npairs,
                                               const size_t          natoms,
                                               const CBasisFunction &a_function,
                                               const CBasisFunction &b_function,
                                               const CBasisFunction &c_function,
                                               const CSimdMatrix    &coordinates,
                                               const CSimdMatrix    &c_coordinates,
                                               CSimdMatrix          &buffer,
                                               const double          omega,
                                               const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_gik_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 424297, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3510 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    simdfunc::compute_pair_exponents(a_function, b_function, coordinates, nmax);

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 424297, 331732, 22470, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

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

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 6, 3, {1, 2, 3,
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                                                            15, 16, 17}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 24, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 69, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 72, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 75, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 78, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 81, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 84, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 87, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 90, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 93, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 96, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 99, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 102, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 105, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 108, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 111, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 114, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 117, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 120, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 123, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 126, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 129, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 132, 0, 3, 39, 40,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 135, 0, 3, 40, 41,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 138, 0, 3, 7, 8,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 144, 0, 3, 8, 9,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 150, 0, 3, 9, 10,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 156, 0, 3, 10, 11,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 162, 0, 3, 11, 12,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 168, 0, 3, 12, 13,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 174, 0, 3, 13, 14,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 180, 0, 3, 14, 15,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 186, 0, 3, 15, 16,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 192, 0, 3, 16, 17,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 198, 0, 3, 17, 18,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 204, 0, 3, 18, 19,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 210, 0, 3, 19, 20,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 216, 0, 3, 20, 21,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 222, 0, 3, 21, 22,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 228, 0, 3, 25, 26,
                                                                       90, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 234, 0, 3, 26, 27,
                                                                       93, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 240, 0, 3, 27, 28,
                                                                       96, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 246, 0, 3, 28, 29,
                                                                       99, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 252, 0, 3, 29, 30,
                                                                       102, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 258, 0, 3, 30, 31,
                                                                       105, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 264, 0, 3, 31, 32,
                                                                       108, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 270, 0, 3, 32, 33,
                                                                       111, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 276, 0, 3, 33, 34,
                                                                       114, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 282, 0, 3, 34, 35,
                                                                       117, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 288, 0, 3, 35, 36,
                                                                       120, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 294, 0, 3, 36, 37,
                                                                       123, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 300, 0, 3, 37, 38,
                                                                       126, 129, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 306, 0, 3, 38, 39,
                                                                       129, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 312, 0, 3, 39, 40,
                                                                       132, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 42, 45,
                                                                       138, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 45, 48,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 48, 51,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 51, 54,
                                                                       156, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 54, 57,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 57, 60,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 60, 63,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 63, 66,
                                                                       180, 186, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 66, 69,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 69, 72,
                                                                       192, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 72, 75,
                                                                       198, 204, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 75, 78,
                                                                       204, 210, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 78, 81,
                                                                       210, 216, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 81, 84,
                                                                       216, 222, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 90, 93,
                                                                       228, 234, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 93, 96,
                                                                       234, 240, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 96, 99,
                                                                       240, 246, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 99,
                                                                       102, 246, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 102,
                                                                       105, 252, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 105,
                                                                       108, 258, 264, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 108,
                                                                       111, 264, 270, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 528, 0, 3, 111,
                                                                       114, 270, 276, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 538, 0, 3, 114,
                                                                       117, 276, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 117,
                                                                       120, 282, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 558, 0, 3, 120,
                                                                       123, 288, 294, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 568, 0, 3, 123,
                                                                       126, 294, 300, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 578, 0, 3, 126,
                                                                       129, 300, 306, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 129,
                                                                       132, 306, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 598, 0, 3, 138,
                                                                       144, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 613, 0, 3, 144,
                                                                       150, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 628, 0, 3, 150,
                                                                       156, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 643, 0, 3, 156,
                                                                       162, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 658, 0, 3, 162,
                                                                       168, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 673, 0, 3, 168,
                                                                       174, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 688, 0, 3, 174,
                                                                       180, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 703, 0, 3, 180,
                                                                       186, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 718, 0, 3, 186,
                                                                       192, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 733, 0, 3, 192,
                                                                       198, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 748, 0, 3, 198,
                                                                       204, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 763, 0, 3, 204,
                                                                       210, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 778, 0, 3, 210,
                                                                       216, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 793, 0, 3, 228,
                                                                       234, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 808, 0, 3, 234,
                                                                       240, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 823, 0, 3, 240,
                                                                       246, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 838, 0, 3, 246,
                                                                       252, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 853, 0, 3, 252,
                                                                       258, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 868, 0, 3, 258,
                                                                       264, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 883, 0, 3, 264,
                                                                       270, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 898, 0, 3, 270,
                                                                       276, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 913, 0, 3, 276,
                                                                       282, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 928, 0, 3, 282,
                                                                       288, 548, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 943, 0, 3, 288,
                                                                       294, 558, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 958, 0, 3, 294,
                                                                       300, 568, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 973, 0, 3, 300,
                                                                       306, 578, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 988, 0, 3, 318,
                                                                       328, 598, 613, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1009, 0, 3, 328,
                                                                       338, 613, 628, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 338,
                                                                       348, 628, 643, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1051, 0, 3, 348,
                                                                       358, 643, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 358,
                                                                       368, 658, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1093, 0, 3, 368,
                                                                       378, 673, 688, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 378,
                                                                       388, 688, 703, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1135, 0, 3, 388,
                                                                       398, 703, 718, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 398,
                                                                       408, 718, 733, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1177, 0, 3, 408,
                                                                       418, 733, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1198, 0, 3, 418,
                                                                       428, 748, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1219, 0, 3, 428,
                                                                       438, 763, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 458,
                                                                       468, 793, 808, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1261, 0, 3, 468,
                                                                       478, 808, 823, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1282, 0, 3, 478,
                                                                       488, 823, 838, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1303, 0, 3, 488,
                                                                       498, 838, 853, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 498,
                                                                       508, 853, 868, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1345, 0, 3, 508,
                                                                       518, 868, 883, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1366, 0, 3, 518,
                                                                       528, 883, 898, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1387, 0, 3, 528,
                                                                       538, 898, 913, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 538,
                                                                       548, 913, 928, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1429, 0, 3, 548,
                                                                       558, 928, 943, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1450, 0, 3, 558,
                                                                       568, 943, 958, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1471, 0, 3, 568,
                                                                       578, 958, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 598,
                                                                       613, 988, 1009, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 613,
                                                                       628, 1009, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 628,
                                                                       643, 1030, 1051, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 643,
                                                                       658, 1051, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 658,
                                                                       673, 1072, 1093, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 673,
                                                                       688, 1093, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 688,
                                                                       703, 1114, 1135, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 703,
                                                                       718, 1135, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 718,
                                                                       733, 1156, 1177, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 733,
                                                                       748, 1177, 1198, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 748,
                                                                       763, 1198, 1219, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 793,
                                                                       808, 1240, 1261, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 808,
                                                                       823, 1261, 1282, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 823,
                                                                       838, 1282, 1303, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 838,
                                                                       853, 1303, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 853,
                                                                       868, 1324, 1345, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 868,
                                                                       883, 1345, 1366, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1968, 0, 3, 883,
                                                                       898, 1366, 1387, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1996, 0, 3, 898,
                                                                       913, 1387, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 913,
                                                                       928, 1408, 1429, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 928,
                                                                       943, 1429, 1450, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 943,
                                                                       958, 1450, 1471, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 988,
                                                                       1009, 1492, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2144, 0, 3, 1009,
                                                                       1030, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2180, 0, 3, 1030,
                                                                       1051, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2216, 0, 3, 1051,
                                                                       1072, 1576, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2252, 0, 3, 1072,
                                                                       1093, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2288, 0, 3, 1093,
                                                                       1114, 1632, 1660, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2324, 0, 3, 1114,
                                                                       1135, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2360, 0, 3, 1135,
                                                                       1156, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2396, 0, 3, 1156,
                                                                       1177, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2432, 0, 3, 1177,
                                                                       1198, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2468, 0, 3, 1240,
                                                                       1261, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2504, 0, 3, 1261,
                                                                       1282, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2540, 0, 3, 1282,
                                                                       1303, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2576, 0, 3, 1303,
                                                                       1324, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2612, 0, 3, 1324,
                                                                       1345, 1912, 1940, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1345,
                                                                       1366, 1940, 1968, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2684, 0, 3, 1366,
                                                                       1387, 1968, 1996, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2720, 0, 3, 1387,
                                                                       1408, 1996, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2756, 0, 3, 1408,
                                                                       1429, 2024, 2052, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2792, 0, 3, 1429,
                                                                       1450, 2052, 2080, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2828, 0, 3, 1492,
                                                                       1520, 2108, 2144, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2873, 0, 3, 1520,
                                                                       1548, 2144, 2180, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2918, 0, 3, 1548,
                                                                       1576, 2180, 2216, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2963, 0, 3, 1576,
                                                                       1604, 2216, 2252, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3008, 0, 3, 1604,
                                                                       1632, 2252, 2288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3053, 0, 3, 1632,
                                                                       1660, 2288, 2324, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3098, 0, 3, 1660,
                                                                       1688, 2324, 2360, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3143, 0, 3, 1688,
                                                                       1716, 2360, 2396, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3188, 0, 3, 1716,
                                                                       1744, 2396, 2432, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3233, 0, 3, 1800,
                                                                       1828, 2468, 2504, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3278, 0, 3, 1828,
                                                                       1856, 2504, 2540, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3323, 0, 3, 1856,
                                                                       1884, 2540, 2576, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3368, 0, 3, 1884,
                                                                       1912, 2576, 2612, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3413, 0, 3, 1912,
                                                                       1940, 2612, 2648, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3458, 0, 3, 1940,
                                                                       1968, 2648, 2684, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3503, 0, 3, 1968,
                                                                       1996, 2684, 2720, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3548, 0, 3, 1996,
                                                                       2024, 2720, 2756, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3593, 0, 3, 2024,
                                                                       2052, 2756, 2792, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3638, 0, 3, 2108,
                                                                       2144, 2828, 2873, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3693, 0, 3, 2144,
                                                                       2180, 2873, 2918, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3748, 0, 3, 2180,
                                                                       2216, 2918, 2963, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3803, 0, 3, 2216,
                                                                       2252, 2963, 3008, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3858, 0, 3, 2252,
                                                                       2288, 3008, 3053, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3913, 0, 3, 2288,
                                                                       2324, 3053, 3098, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2324,
                                                                       2360, 3098, 3143, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4023, 0, 3, 2360,
                                                                       2396, 3143, 3188, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4078, 0, 3, 2468,
                                                                       2504, 3233, 3278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4133, 0, 3, 2504,
                                                                       2540, 3278, 3323, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4188, 0, 3, 2540,
                                                                       2576, 3323, 3368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4243, 0, 3, 2576,
                                                                       2612, 3368, 3413, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 2612,
                                                                       2648, 3413, 3458, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4353, 0, 3, 2648,
                                                                       2684, 3458, 3503, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4408, 0, 3, 2684,
                                                                       2720, 3503, 3548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4463, 0, 3, 2720,
                                                                       2756, 3548, 3593, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4518, 0, 3, 2828,
                                                                       2873, 3638, 3693, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4584, 0, 3, 2873,
                                                                       2918, 3693, 3748, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4650, 0, 3, 2918,
                                                                       2963, 3748, 3803, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4716, 0, 3, 2963,
                                                                       3008, 3803, 3858, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4782, 0, 3, 3008,
                                                                       3053, 3858, 3913, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4848, 0, 3, 3053,
                                                                       3098, 3913, 3968, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4914, 0, 3, 3098,
                                                                       3143, 3968, 4023, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4980, 0, 3, 3233,
                                                                       3278, 4078, 4133, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5046, 0, 3, 3278,
                                                                       3323, 4133, 4188, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5112, 0, 3, 3323,
                                                                       3368, 4188, 4243, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5178, 0, 3, 3368,
                                                                       3413, 4243, 4298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5244, 0, 3, 3413,
                                                                       3458, 4298, 4353, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5310, 0, 3, 3458,
                                                                       3503, 4353, 4408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5376, 0, 3, 3503,
                                                                       3548, 4408, 4463, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5442, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5445, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5448, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5451, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5454, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5457, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5460, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5463, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5466, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5469, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5472, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5475, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5478, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5481, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5484, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5487, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5490, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5493, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5496, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5499, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5502, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5505, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5508, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5511, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5514, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5517, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5520, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5523, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5526, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5529, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5532, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5535, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5538, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5541, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5544, 3, 9, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5553, 3, 10, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5562, 3, 11, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5571, 3, 12, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5580, 3, 13, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5589, 3, 14, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5598, 3, 15, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5607, 3, 16, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5616, 3, 17, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5625, 3, 18, 75,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5634, 3, 19, 78,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5643, 3, 20, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5652, 3, 21, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5661, 3, 22, 87,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5670, 3, 27, 96,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5679, 3, 28, 99,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5688, 3, 29, 102,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5697, 3, 30, 105,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5706, 3, 31, 108,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5715, 3, 32, 111,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5724, 3, 33, 114,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5733, 3, 34, 117,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5742, 3, 35, 120,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5751, 3, 36, 123,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5760, 3, 37, 126,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5769, 3, 38, 129,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5778, 3, 39, 132,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5787, 3, 40, 135,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5796, 3, 42, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5814, 3, 45, 144,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5832, 3, 48, 150,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5850, 3, 51, 156,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5868, 3, 54, 162,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5886, 3, 57, 168,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5904, 3, 60, 174,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5922, 3, 63, 180,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5940, 3, 66, 186,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5958, 3, 69, 192,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5976, 3, 72, 198,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5994, 3, 75, 204,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6012, 3, 78, 210,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6030, 3, 81, 216,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6048, 3, 84, 222,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6066, 3, 90, 228,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6084, 3, 93, 234,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6102, 3, 96, 240,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6120, 3, 99, 246,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6138, 3, 102, 252,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6156, 3, 105, 258,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6174, 3, 108, 264,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6192, 3, 111, 270,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6210, 3, 114, 276,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6228, 3, 117, 282,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6246, 3, 120, 288,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6264, 3, 123, 294,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6282, 3, 126, 300,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6300, 3, 129, 306,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6318, 3, 132, 312,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6336, 3, 138, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6366, 3, 144, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6396, 3, 150, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6426, 3, 156, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6456, 3, 162, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6486, 3, 168, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6516, 3, 174, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6546, 3, 180, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6576, 3, 186, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6606, 3, 192, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6636, 3, 198, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6666, 3, 204, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6696, 3, 210, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6726, 3, 216, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6756, 3, 228, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6786, 3, 234, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6816, 3, 240, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6846, 3, 246, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6876, 3, 252, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6906, 3, 258, 508,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6936, 3, 264, 518,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6966, 3, 270, 528,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6996, 3, 276, 538,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7026, 3, 282, 548,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7056, 3, 288, 558,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7086, 3, 294, 568,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7116, 3, 300, 578,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7146, 3, 306, 588,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7176, 3, 318, 598,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7221, 3, 328, 613,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7266, 3, 338, 628,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7311, 3, 348, 643,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7356, 3, 358, 658,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7401, 3, 368, 673,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7446, 3, 378, 688,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7491, 3, 388, 703,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7536, 3, 398, 718,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7581, 3, 408, 733,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7626, 3, 418, 748,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7671, 3, 428, 763,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7716, 3, 438, 778,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7761, 3, 458, 793,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7806, 3, 468, 808,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7851, 3, 478, 823,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7896, 3, 488, 838,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7941, 3, 498, 853,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7986, 3, 508, 868,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8031, 3, 518, 883,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8076, 3, 528, 898,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8121, 3, 538, 913,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8166, 3, 548, 928,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8211, 3, 558, 943,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8256, 3, 568, 958,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8301, 3, 578, 973,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8346, 3, 598, 988,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8409, 3, 613,
                                                                       1009, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8472, 3, 628,
                                                                       1030, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8535, 3, 643,
                                                                       1051, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8598, 3, 658,
                                                                       1072, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8661, 3, 673,
                                                                       1093, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8724, 3, 688,
                                                                       1114, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8787, 3, 703,
                                                                       1135, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8850, 3, 718,
                                                                       1156, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8913, 3, 733,
                                                                       1177, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8976, 3, 748,
                                                                       1198, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9039, 3, 763,
                                                                       1219, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9102, 3, 793,
                                                                       1240, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9165, 3, 808,
                                                                       1261, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9228, 3, 823,
                                                                       1282, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9291, 3, 838,
                                                                       1303, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9354, 3, 853,
                                                                       1324, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9417, 3, 868,
                                                                       1345, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9480, 3, 883,
                                                                       1366, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9543, 3, 898,
                                                                       1387, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9606, 3, 913,
                                                                       1408, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9669, 3, 928,
                                                                       1429, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9732, 3, 943,
                                                                       1450, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9795, 3, 958,
                                                                       1471, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9858, 3, 988,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9942, 3, 1009,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10026, 3, 1030,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10110, 3, 1051,
                                                                       1576, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10194, 3, 1072,
                                                                       1604, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10278, 3, 1093,
                                                                       1632, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10362, 3, 1114,
                                                                       1660, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10446, 3, 1135,
                                                                       1688, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10530, 3, 1156,
                                                                       1716, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10614, 3, 1177,
                                                                       1744, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10698, 3, 1198,
                                                                       1772, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10782, 3, 1240,
                                                                       1800, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10866, 3, 1261,
                                                                       1828, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10950, 3, 1282,
                                                                       1856, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11034, 3, 1303,
                                                                       1884, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11118, 3, 1324,
                                                                       1912, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11202, 3, 1345,
                                                                       1940, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11286, 3, 1366,
                                                                       1968, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11370, 3, 1387,
                                                                       1996, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11454, 3, 1408,
                                                                       2024, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11538, 3, 1429,
                                                                       2052, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11622, 3, 1450,
                                                                       2080, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11706, 3, 1492,
                                                                       2108, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11814, 3, 1520,
                                                                       2144, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11922, 3, 1548,
                                                                       2180, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12030, 3, 1576,
                                                                       2216, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12138, 3, 1604,
                                                                       2252, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12246, 3, 1632,
                                                                       2288, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12354, 3, 1660,
                                                                       2324, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12462, 3, 1688,
                                                                       2360, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12570, 3, 1716,
                                                                       2396, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12678, 3, 1744,
                                                                       2432, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12786, 3, 1800,
                                                                       2468, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12894, 3, 1828,
                                                                       2504, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13002, 3, 1856,
                                                                       2540, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13110, 3, 1884,
                                                                       2576, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13218, 3, 1912,
                                                                       2612, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13326, 3, 1940,
                                                                       2648, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13434, 3, 1968,
                                                                       2684, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13542, 3, 1996,
                                                                       2720, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13650, 3, 2024,
                                                                       2756, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13758, 3, 2052,
                                                                       2792, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13866, 3, 2108,
                                                                       2828, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14001, 3, 2144,
                                                                       2873, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14136, 3, 2180,
                                                                       2918, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14271, 3, 2216,
                                                                       2963, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14406, 3, 2252,
                                                                       3008, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14541, 3, 2288,
                                                                       3053, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14676, 3, 2324,
                                                                       3098, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14811, 3, 2360,
                                                                       3143, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14946, 3, 2396,
                                                                       3188, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15081, 3, 2468,
                                                                       3233, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15216, 3, 2504,
                                                                       3278, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15351, 3, 2540,
                                                                       3323, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15486, 3, 2576,
                                                                       3368, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15621, 3, 2612,
                                                                       3413, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15756, 3, 2648,
                                                                       3458, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15891, 3, 2684,
                                                                       3503, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16026, 3, 2720,
                                                                       3548, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16161, 3, 2756,
                                                                       3593, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16296, 3, 2828,
                                                                       3638, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16461, 3, 2873,
                                                                       3693, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16626, 3, 2918,
                                                                       3748, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16791, 3, 2963,
                                                                       3803, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16956, 3, 3008,
                                                                       3858, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17121, 3, 3053,
                                                                       3913, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17286, 3, 3098,
                                                                       3968, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17451, 3, 3143,
                                                                       4023, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17616, 3, 3233,
                                                                       4078, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17781, 3, 3278,
                                                                       4133, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17946, 3, 3323,
                                                                       4188, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18111, 3, 3368,
                                                                       4243, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18276, 3, 3413,
                                                                       4298, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18441, 3, 3458,
                                                                       4353, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18606, 3, 3503,
                                                                       4408, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18771, 3, 3548,
                                                                       4463, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 18936, 3, 3638,
                                                                       4518, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 19134, 3, 3693,
                                                                       4584, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 19332, 3, 3748,
                                                                       4650, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 19530, 3, 3803,
                                                                       4716, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 19728, 3, 3858,
                                                                       4782, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 19926, 3, 3913,
                                                                       4848, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20124, 3, 3968,
                                                                       4914, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20322, 3, 4078,
                                                                       4980, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20520, 3, 4133,
                                                                       5046, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20718, 3, 4188,
                                                                       5112, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20916, 3, 4243,
                                                                       5178, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 21114, 3, 4298,
                                                                       5244, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 21312, 3, 4353,
                                                                       5310, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 21510, 3, 4408,
                                                                       5376, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21708, 3, 7, 8,
                                                                       5448, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21714, 3, 8, 9,
                                                                       5451, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21720, 3, 9, 10,
                                                                       5454, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21726, 3, 10, 11,
                                                                       5457, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21732, 3, 11, 12,
                                                                       5460, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21738, 3, 12, 13,
                                                                       5463, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21744, 3, 13, 14,
                                                                       5466, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21750, 3, 14, 15,
                                                                       5469, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21756, 3, 15, 16,
                                                                       5472, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21762, 3, 16, 17,
                                                                       5475, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21768, 3, 17, 18,
                                                                       5478, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21774, 3, 18, 19,
                                                                       5481, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21780, 3, 19, 20,
                                                                       5484, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21786, 3, 20, 21,
                                                                       5487, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21792, 3, 21, 22,
                                                                       5490, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21798, 3, 25, 26,
                                                                       5499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21804, 3, 26, 27,
                                                                       5502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21810, 3, 27, 28,
                                                                       5505, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21816, 3, 28, 29,
                                                                       5508, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21822, 3, 29, 30,
                                                                       5511, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21828, 3, 30, 31,
                                                                       5514, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21834, 3, 31, 32,
                                                                       5517, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21840, 3, 32, 33,
                                                                       5520, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21846, 3, 33, 34,
                                                                       5523, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21852, 3, 34, 35,
                                                                       5526, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21858, 3, 35, 36,
                                                                       5529, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21864, 3, 36, 37,
                                                                       5532, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21870, 3, 37, 38,
                                                                       5535, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21876, 3, 38, 39,
                                                                       5538, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 21882, 3, 39, 40,
                                                                       5541, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 21888, 0, 3,
                                                                       21708, 5448, 21714, 5544,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 21906, 0, 3,
                                                                       21714, 5451, 21720, 5553,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 21924, 0, 3,
                                                                       21720, 5454, 21726, 5562,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 21942, 0, 3,
                                                                       21726, 5457, 21732, 5571,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 21960, 0, 3,
                                                                       21732, 5460, 21738, 5580,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 21978, 0, 3,
                                                                       21738, 5463, 21744, 5589,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 21996, 0, 3,
                                                                       21744, 5466, 21750, 5598,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22014, 0, 3,
                                                                       21750, 5469, 21756, 5607,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22032, 0, 3,
                                                                       21756, 5472, 21762, 5616,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22050, 0, 3,
                                                                       21762, 5475, 21768, 5625,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22068, 0, 3,
                                                                       21768, 5478, 21774, 5634,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22086, 0, 3,
                                                                       21774, 5481, 21780, 5643,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22104, 0, 3,
                                                                       21780, 5484, 21786, 5652,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22122, 0, 3,
                                                                       21786, 5487, 21792, 5661,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22140, 0, 3,
                                                                       21798, 5499, 21804, 5670,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22158, 0, 3,
                                                                       21804, 5502, 21810, 5679,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22176, 0, 3,
                                                                       21810, 5505, 21816, 5688,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22194, 0, 3,
                                                                       21816, 5508, 21822, 5697,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22212, 0, 3,
                                                                       21822, 5511, 21828, 5706,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22230, 0, 3,
                                                                       21828, 5514, 21834, 5715,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22248, 0, 3,
                                                                       21834, 5517, 21840, 5724,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22266, 0, 3,
                                                                       21840, 5520, 21846, 5733,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22284, 0, 3,
                                                                       21846, 5523, 21852, 5742,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22302, 0, 3,
                                                                       21852, 5526, 21858, 5751,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22320, 0, 3,
                                                                       21858, 5529, 21864, 5760,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22338, 0, 3,
                                                                       21864, 5532, 21870, 5769,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22356, 0, 3,
                                                                       21870, 5535, 21876, 5778,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22374, 0, 3,
                                                                       21876, 5538, 21882, 5787,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22392, 0, 3,
                                                                       21888, 5544, 21906, 138,
                                                                       144, 5832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22428, 0, 3,
                                                                       21906, 5553, 21924, 144,
                                                                       150, 5850, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22464, 0, 3,
                                                                       21924, 5562, 21942, 150,
                                                                       156, 5868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22500, 0, 3,
                                                                       21942, 5571, 21960, 156,
                                                                       162, 5886, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22536, 0, 3,
                                                                       21960, 5580, 21978, 162,
                                                                       168, 5904, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22572, 0, 3,
                                                                       21978, 5589, 21996, 168,
                                                                       174, 5922, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22608, 0, 3,
                                                                       21996, 5598, 22014, 174,
                                                                       180, 5940, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22644, 0, 3,
                                                                       22014, 5607, 22032, 180,
                                                                       186, 5958, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22680, 0, 3,
                                                                       22032, 5616, 22050, 186,
                                                                       192, 5976, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22716, 0, 3,
                                                                       22050, 5625, 22068, 192,
                                                                       198, 5994, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22752, 0, 3,
                                                                       22068, 5634, 22086, 198,
                                                                       204, 6012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22788, 0, 3,
                                                                       22086, 5643, 22104, 204,
                                                                       210, 6030, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22824, 0, 3,
                                                                       22104, 5652, 22122, 210,
                                                                       216, 6048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22860, 0, 3,
                                                                       22140, 5670, 22158, 228,
                                                                       234, 6102, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22896, 0, 3,
                                                                       22158, 5679, 22176, 234,
                                                                       240, 6120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22932, 0, 3,
                                                                       22176, 5688, 22194, 240,
                                                                       246, 6138, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 22968, 0, 3,
                                                                       22194, 5697, 22212, 246,
                                                                       252, 6156, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23004, 0, 3,
                                                                       22212, 5706, 22230, 252,
                                                                       258, 6174, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23040, 0, 3,
                                                                       22230, 5715, 22248, 258,
                                                                       264, 6192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23076, 0, 3,
                                                                       22248, 5724, 22266, 264,
                                                                       270, 6210, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23112, 0, 3,
                                                                       22266, 5733, 22284, 270,
                                                                       276, 6228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23148, 0, 3,
                                                                       22284, 5742, 22302, 276,
                                                                       282, 6246, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23184, 0, 3,
                                                                       22302, 5751, 22320, 282,
                                                                       288, 6264, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23220, 0, 3,
                                                                       22320, 5760, 22338, 288,
                                                                       294, 6282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23256, 0, 3,
                                                                       22338, 5769, 22356, 294,
                                                                       300, 6300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23292, 0, 3,
                                                                       22356, 5778, 22374, 300,
                                                                       306, 6318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23328, 0, 3,
                                                                       22392, 5832, 22428, 318,
                                                                       328, 6396, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23388, 0, 3,
                                                                       22428, 5850, 22464, 328,
                                                                       338, 6426, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23448, 0, 3,
                                                                       22464, 5868, 22500, 338,
                                                                       348, 6456, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23508, 0, 3,
                                                                       22500, 5886, 22536, 348,
                                                                       358, 6486, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23568, 0, 3,
                                                                       22536, 5904, 22572, 358,
                                                                       368, 6516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23628, 0, 3,
                                                                       22572, 5922, 22608, 368,
                                                                       378, 6546, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23688, 0, 3,
                                                                       22608, 5940, 22644, 378,
                                                                       388, 6576, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23748, 0, 3,
                                                                       22644, 5958, 22680, 388,
                                                                       398, 6606, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23808, 0, 3,
                                                                       22680, 5976, 22716, 398,
                                                                       408, 6636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23868, 0, 3,
                                                                       22716, 5994, 22752, 408,
                                                                       418, 6666, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23928, 0, 3,
                                                                       22752, 6012, 22788, 418,
                                                                       428, 6696, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 23988, 0, 3,
                                                                       22788, 6030, 22824, 428,
                                                                       438, 6726, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24048, 0, 3,
                                                                       22860, 6102, 22896, 458,
                                                                       468, 6816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24108, 0, 3,
                                                                       22896, 6120, 22932, 468,
                                                                       478, 6846, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24168, 0, 3,
                                                                       22932, 6138, 22968, 478,
                                                                       488, 6876, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24228, 0, 3,
                                                                       22968, 6156, 23004, 488,
                                                                       498, 6906, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24288, 0, 3,
                                                                       23004, 6174, 23040, 498,
                                                                       508, 6936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24348, 0, 3,
                                                                       23040, 6192, 23076, 508,
                                                                       518, 6966, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24408, 0, 3,
                                                                       23076, 6210, 23112, 518,
                                                                       528, 6996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24468, 0, 3,
                                                                       23112, 6228, 23148, 528,
                                                                       538, 7026, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24528, 0, 3,
                                                                       23148, 6246, 23184, 538,
                                                                       548, 7056, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24588, 0, 3,
                                                                       23184, 6264, 23220, 548,
                                                                       558, 7086, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24648, 0, 3,
                                                                       23220, 6282, 23256, 558,
                                                                       568, 7116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24708, 0, 3,
                                                                       23256, 6300, 23292, 568,
                                                                       578, 7146, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 24768, 0, 3,
                                                                       23328, 6396, 23388, 598,
                                                                       613, 7266, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 24858, 0, 3,
                                                                       23388, 6426, 23448, 613,
                                                                       628, 7311, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 24948, 0, 3,
                                                                       23448, 6456, 23508, 628,
                                                                       643, 7356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25038, 0, 3,
                                                                       23508, 6486, 23568, 643,
                                                                       658, 7401, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25128, 0, 3,
                                                                       23568, 6516, 23628, 658,
                                                                       673, 7446, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25218, 0, 3,
                                                                       23628, 6546, 23688, 673,
                                                                       688, 7491, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25308, 0, 3,
                                                                       23688, 6576, 23748, 688,
                                                                       703, 7536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25398, 0, 3,
                                                                       23748, 6606, 23808, 703,
                                                                       718, 7581, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25488, 0, 3,
                                                                       23808, 6636, 23868, 718,
                                                                       733, 7626, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25578, 0, 3,
                                                                       23868, 6666, 23928, 733,
                                                                       748, 7671, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25668, 0, 3,
                                                                       23928, 6696, 23988, 748,
                                                                       763, 7716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25758, 0, 3,
                                                                       24048, 6816, 24108, 793,
                                                                       808, 7851, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25848, 0, 3,
                                                                       24108, 6846, 24168, 808,
                                                                       823, 7896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25938, 0, 3,
                                                                       24168, 6876, 24228, 823,
                                                                       838, 7941, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26028, 0, 3,
                                                                       24228, 6906, 24288, 838,
                                                                       853, 7986, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26118, 0, 3,
                                                                       24288, 6936, 24348, 853,
                                                                       868, 8031, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26208, 0, 3,
                                                                       24348, 6966, 24408, 868,
                                                                       883, 8076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26298, 0, 3,
                                                                       24408, 6996, 24468, 883,
                                                                       898, 8121, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26388, 0, 3,
                                                                       24468, 7026, 24528, 898,
                                                                       913, 8166, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26478, 0, 3,
                                                                       24528, 7056, 24588, 913,
                                                                       928, 8211, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26568, 0, 3,
                                                                       24588, 7086, 24648, 928,
                                                                       943, 8256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26658, 0, 3,
                                                                       24648, 7116, 24708, 943,
                                                                       958, 8301, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 26748, 0, 3,
                                                                       24768, 7266, 24858, 988,
                                                                       1009, 8472, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 26874, 0, 3,
                                                                       24858, 7311, 24948, 1009,
                                                                       1030, 8535, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27000, 0, 3,
                                                                       24948, 7356, 25038, 1030,
                                                                       1051, 8598, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27126, 0, 3,
                                                                       25038, 7401, 25128, 1051,
                                                                       1072, 8661, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27252, 0, 3,
                                                                       25128, 7446, 25218, 1072,
                                                                       1093, 8724, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27378, 0, 3,
                                                                       25218, 7491, 25308, 1093,
                                                                       1114, 8787, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27504, 0, 3,
                                                                       25308, 7536, 25398, 1114,
                                                                       1135, 8850, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27630, 0, 3,
                                                                       25398, 7581, 25488, 1135,
                                                                       1156, 8913, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27756, 0, 3,
                                                                       25488, 7626, 25578, 1156,
                                                                       1177, 8976, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27882, 0, 3,
                                                                       25578, 7671, 25668, 1177,
                                                                       1198, 9039, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28008, 0, 3,
                                                                       25758, 7851, 25848, 1240,
                                                                       1261, 9228, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28134, 0, 3,
                                                                       25848, 7896, 25938, 1261,
                                                                       1282, 9291, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28260, 0, 3,
                                                                       25938, 7941, 26028, 1282,
                                                                       1303, 9354, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28386, 0, 3,
                                                                       26028, 7986, 26118, 1303,
                                                                       1324, 9417, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28512, 0, 3,
                                                                       26118, 8031, 26208, 1324,
                                                                       1345, 9480, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28638, 0, 3,
                                                                       26208, 8076, 26298, 1345,
                                                                       1366, 9543, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28764, 0, 3,
                                                                       26298, 8121, 26388, 1366,
                                                                       1387, 9606, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28890, 0, 3,
                                                                       26388, 8166, 26478, 1387,
                                                                       1408, 9669, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29016, 0, 3,
                                                                       26478, 8211, 26568, 1408,
                                                                       1429, 9732, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29142, 0, 3,
                                                                       26568, 8256, 26658, 1429,
                                                                       1450, 9795, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 29268, 0, 3,
                                                                       26748, 8472, 26874, 1492,
                                                                       1520, 10026, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 29436, 0, 3,
                                                                       26874, 8535, 27000, 1520,
                                                                       1548, 10110, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 29604, 0, 3,
                                                                       27000, 8598, 27126, 1548,
                                                                       1576, 10194, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 29772, 0, 3,
                                                                       27126, 8661, 27252, 1576,
                                                                       1604, 10278, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 29940, 0, 3,
                                                                       27252, 8724, 27378, 1604,
                                                                       1632, 10362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30108, 0, 3,
                                                                       27378, 8787, 27504, 1632,
                                                                       1660, 10446, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30276, 0, 3,
                                                                       27504, 8850, 27630, 1660,
                                                                       1688, 10530, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30444, 0, 3,
                                                                       27630, 8913, 27756, 1688,
                                                                       1716, 10614, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30612, 0, 3,
                                                                       27756, 8976, 27882, 1716,
                                                                       1744, 10698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30780, 0, 3,
                                                                       28008, 9228, 28134, 1800,
                                                                       1828, 10950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30948, 0, 3,
                                                                       28134, 9291, 28260, 1828,
                                                                       1856, 11034, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31116, 0, 3,
                                                                       28260, 9354, 28386, 1856,
                                                                       1884, 11118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31284, 0, 3,
                                                                       28386, 9417, 28512, 1884,
                                                                       1912, 11202, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31452, 0, 3,
                                                                       28512, 9480, 28638, 1912,
                                                                       1940, 11286, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31620, 0, 3,
                                                                       28638, 9543, 28764, 1940,
                                                                       1968, 11370, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31788, 0, 3,
                                                                       28764, 9606, 28890, 1968,
                                                                       1996, 11454, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31956, 0, 3,
                                                                       28890, 9669, 29016, 1996,
                                                                       2024, 11538, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 32124, 0, 3,
                                                                       29016, 9732, 29142, 2024,
                                                                       2052, 11622, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 32292, 0, 3,
                                                                       29268, 10026, 29436, 2108,
                                                                       2144, 11922, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 32508, 0, 3,
                                                                       29436, 10110, 29604, 2144,
                                                                       2180, 12030, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 32724, 0, 3,
                                                                       29604, 10194, 29772, 2180,
                                                                       2216, 12138, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 32940, 0, 3,
                                                                       29772, 10278, 29940, 2216,
                                                                       2252, 12246, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 33156, 0, 3,
                                                                       29940, 10362, 30108, 2252,
                                                                       2288, 12354, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 33372, 0, 3,
                                                                       30108, 10446, 30276, 2288,
                                                                       2324, 12462, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 33588, 0, 3,
                                                                       30276, 10530, 30444, 2324,
                                                                       2360, 12570, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 33804, 0, 3,
                                                                       30444, 10614, 30612, 2360,
                                                                       2396, 12678, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 34020, 0, 3,
                                                                       30780, 10950, 30948, 2468,
                                                                       2504, 13002, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 34236, 0, 3,
                                                                       30948, 11034, 31116, 2504,
                                                                       2540, 13110, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 34452, 0, 3,
                                                                       31116, 11118, 31284, 2540,
                                                                       2576, 13218, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 34668, 0, 3,
                                                                       31284, 11202, 31452, 2576,
                                                                       2612, 13326, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 34884, 0, 3,
                                                                       31452, 11286, 31620, 2612,
                                                                       2648, 13434, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 35100, 0, 3,
                                                                       31620, 11370, 31788, 2648,
                                                                       2684, 13542, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 35316, 0, 3,
                                                                       31788, 11454, 31956, 2684,
                                                                       2720, 13650, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 35532, 0, 3,
                                                                       31956, 11538, 32124, 2720,
                                                                       2756, 13758, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 35748, 0, 3,
                                                                       32292, 11922, 32508, 2828,
                                                                       2873, 14136, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 36018, 0, 3,
                                                                       32508, 12030, 32724, 2873,
                                                                       2918, 14271, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 36288, 0, 3,
                                                                       32724, 12138, 32940, 2918,
                                                                       2963, 14406, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 36558, 0, 3,
                                                                       32940, 12246, 33156, 2963,
                                                                       3008, 14541, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 36828, 0, 3,
                                                                       33156, 12354, 33372, 3008,
                                                                       3053, 14676, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 37098, 0, 3,
                                                                       33372, 12462, 33588, 3053,
                                                                       3098, 14811, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 37368, 0, 3,
                                                                       33588, 12570, 33804, 3098,
                                                                       3143, 14946, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 37638, 0, 3,
                                                                       34020, 13002, 34236, 3233,
                                                                       3278, 15351, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 37908, 0, 3,
                                                                       34236, 13110, 34452, 3278,
                                                                       3323, 15486, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 38178, 0, 3,
                                                                       34452, 13218, 34668, 3323,
                                                                       3368, 15621, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 38448, 0, 3,
                                                                       34668, 13326, 34884, 3368,
                                                                       3413, 15756, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 38718, 0, 3,
                                                                       34884, 13434, 35100, 3413,
                                                                       3458, 15891, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 38988, 0, 3,
                                                                       35100, 13542, 35316, 3458,
                                                                       3503, 16026, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 39258, 0, 3,
                                                                       35316, 13650, 35532, 3503,
                                                                       3548, 16161, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 39528, 0, 3,
                                                                       35748, 14136, 36018, 3638,
                                                                       3693, 16626, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 39858, 0, 3,
                                                                       36018, 14271, 36288, 3693,
                                                                       3748, 16791, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 40188, 0, 3,
                                                                       36288, 14406, 36558, 3748,
                                                                       3803, 16956, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 40518, 0, 3,
                                                                       36558, 14541, 36828, 3803,
                                                                       3858, 17121, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 40848, 0, 3,
                                                                       36828, 14676, 37098, 3858,
                                                                       3913, 17286, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 41178, 0, 3,
                                                                       37098, 14811, 37368, 3913,
                                                                       3968, 17451, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 41508, 0, 3,
                                                                       37638, 15351, 37908, 4078,
                                                                       4133, 17946, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 41838, 0, 3,
                                                                       37908, 15486, 38178, 4133,
                                                                       4188, 18111, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 42168, 0, 3,
                                                                       38178, 15621, 38448, 4188,
                                                                       4243, 18276, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 42498, 0, 3,
                                                                       38448, 15756, 38718, 4243,
                                                                       4298, 18441, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 42828, 0, 3,
                                                                       38718, 15891, 38988, 4298,
                                                                       4353, 18606, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 43158, 0, 3,
                                                                       38988, 16026, 39258, 4353,
                                                                       4408, 18771, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 43488, 0, 3,
                                                                       39528, 16626, 39858, 4518,
                                                                       4584, 19332, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 43884, 0, 3,
                                                                       39858, 16791, 40188, 4584,
                                                                       4650, 19530, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 44280, 0, 3,
                                                                       40188, 16956, 40518, 4650,
                                                                       4716, 19728, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 44676, 0, 3,
                                                                       40518, 17121, 40848, 4716,
                                                                       4782, 19926, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 45072, 0, 3,
                                                                       40848, 17286, 41178, 4782,
                                                                       4848, 20124, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 45468, 0, 3,
                                                                       41508, 17946, 41838, 4980,
                                                                       5046, 20718, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 45864, 0, 3,
                                                                       41838, 18111, 42168, 5046,
                                                                       5112, 20916, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 46260, 0, 3,
                                                                       42168, 18276, 42498, 5112,
                                                                       5178, 21114, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 46656, 0, 3,
                                                                       42498, 18441, 42828, 5178,
                                                                       5244, 21312, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 47052, 0, 3,
                                                                       42828, 18606, 43158, 5244,
                                                                       5310, 21510, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47448, 3, 5442,
                                                                       5445, 21708, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47458, 3, 5445,
                                                                       5448, 21714, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47468, 3, 5448,
                                                                       5451, 21720, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47478, 3, 5451,
                                                                       5454, 21726, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47488, 3, 5454,
                                                                       5457, 21732, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47498, 3, 5457,
                                                                       5460, 21738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47508, 3, 5460,
                                                                       5463, 21744, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47518, 3, 5463,
                                                                       5466, 21750, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47528, 3, 5466,
                                                                       5469, 21756, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47538, 3, 5469,
                                                                       5472, 21762, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47548, 3, 5472,
                                                                       5475, 21768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47558, 3, 5475,
                                                                       5478, 21774, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47568, 3, 5478,
                                                                       5481, 21780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47578, 3, 5481,
                                                                       5484, 21786, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47588, 3, 5484,
                                                                       5487, 21792, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47598, 3, 5493,
                                                                       5496, 21798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47608, 3, 5496,
                                                                       5499, 21804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47618, 3, 5499,
                                                                       5502, 21810, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47628, 3, 5502,
                                                                       5505, 21816, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47638, 3, 5505,
                                                                       5508, 21822, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47648, 3, 5508,
                                                                       5511, 21828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47658, 3, 5511,
                                                                       5514, 21834, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47668, 3, 5514,
                                                                       5517, 21840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47678, 3, 5517,
                                                                       5520, 21846, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47688, 3, 5520,
                                                                       5523, 21852, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47698, 3, 5523,
                                                                       5526, 21858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47708, 3, 5526,
                                                                       5529, 21864, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47718, 3, 5529,
                                                                       5532, 21870, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47728, 3, 5532,
                                                                       5535, 21876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47738, 3, 5535,
                                                                       5538, 21882, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47748, 0, 3,
                                                                       47448, 21708, 47458,
                                                                       21888, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47778, 0, 3,
                                                                       47458, 21714, 47468,
                                                                       21906, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47808, 0, 3,
                                                                       47468, 21720, 47478,
                                                                       21924, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47838, 0, 3,
                                                                       47478, 21726, 47488,
                                                                       21942, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47868, 0, 3,
                                                                       47488, 21732, 47498,
                                                                       21960, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47898, 0, 3,
                                                                       47498, 21738, 47508,
                                                                       21978, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47928, 0, 3,
                                                                       47508, 21744, 47518,
                                                                       21996, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47958, 0, 3,
                                                                       47518, 21750, 47528,
                                                                       22014, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47988, 0, 3,
                                                                       47528, 21756, 47538,
                                                                       22032, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48018, 0, 3,
                                                                       47538, 21762, 47548,
                                                                       22050, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48048, 0, 3,
                                                                       47548, 21768, 47558,
                                                                       22068, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48078, 0, 3,
                                                                       47558, 21774, 47568,
                                                                       22086, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48108, 0, 3,
                                                                       47568, 21780, 47578,
                                                                       22104, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48138, 0, 3,
                                                                       47578, 21786, 47588,
                                                                       22122, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48168, 0, 3,
                                                                       47598, 21798, 47608,
                                                                       22140, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48198, 0, 3,
                                                                       47608, 21804, 47618,
                                                                       22158, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48228, 0, 3,
                                                                       47618, 21810, 47628,
                                                                       22176, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48258, 0, 3,
                                                                       47628, 21816, 47638,
                                                                       22194, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48288, 0, 3,
                                                                       47638, 21822, 47648,
                                                                       22212, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48318, 0, 3,
                                                                       47648, 21828, 47658,
                                                                       22230, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48348, 0, 3,
                                                                       47658, 21834, 47668,
                                                                       22248, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48378, 0, 3,
                                                                       47668, 21840, 47678,
                                                                       22266, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48408, 0, 3,
                                                                       47678, 21846, 47688,
                                                                       22284, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48438, 0, 3,
                                                                       47688, 21852, 47698,
                                                                       22302, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48468, 0, 3,
                                                                       47698, 21858, 47708,
                                                                       22320, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48498, 0, 3,
                                                                       47708, 21864, 47718,
                                                                       22338, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48528, 0, 3,
                                                                       47718, 21870, 47728,
                                                                       22356, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48558, 0, 3,
                                                                       47728, 21876, 47738,
                                                                       22374, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48588, 0, 3,
                                                                       47748, 21888, 47778, 5796,
                                                                       5814, 22392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48648, 0, 3,
                                                                       47778, 21906, 47808, 5814,
                                                                       5832, 22428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48708, 0, 3,
                                                                       47808, 21924, 47838, 5832,
                                                                       5850, 22464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48768, 0, 3,
                                                                       47838, 21942, 47868, 5850,
                                                                       5868, 22500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48828, 0, 3,
                                                                       47868, 21960, 47898, 5868,
                                                                       5886, 22536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48888, 0, 3,
                                                                       47898, 21978, 47928, 5886,
                                                                       5904, 22572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48948, 0, 3,
                                                                       47928, 21996, 47958, 5904,
                                                                       5922, 22608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49008, 0, 3,
                                                                       47958, 22014, 47988, 5922,
                                                                       5940, 22644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49068, 0, 3,
                                                                       47988, 22032, 48018, 5940,
                                                                       5958, 22680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49128, 0, 3,
                                                                       48018, 22050, 48048, 5958,
                                                                       5976, 22716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49188, 0, 3,
                                                                       48048, 22068, 48078, 5976,
                                                                       5994, 22752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49248, 0, 3,
                                                                       48078, 22086, 48108, 5994,
                                                                       6012, 22788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49308, 0, 3,
                                                                       48108, 22104, 48138, 6012,
                                                                       6030, 22824, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49368, 0, 3,
                                                                       48168, 22140, 48198, 6066,
                                                                       6084, 22860, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49428, 0, 3,
                                                                       48198, 22158, 48228, 6084,
                                                                       6102, 22896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49488, 0, 3,
                                                                       48228, 22176, 48258, 6102,
                                                                       6120, 22932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49548, 0, 3,
                                                                       48258, 22194, 48288, 6120,
                                                                       6138, 22968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49608, 0, 3,
                                                                       48288, 22212, 48318, 6138,
                                                                       6156, 23004, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49668, 0, 3,
                                                                       48318, 22230, 48348, 6156,
                                                                       6174, 23040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49728, 0, 3,
                                                                       48348, 22248, 48378, 6174,
                                                                       6192, 23076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49788, 0, 3,
                                                                       48378, 22266, 48408, 6192,
                                                                       6210, 23112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49848, 0, 3,
                                                                       48408, 22284, 48438, 6210,
                                                                       6228, 23148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49908, 0, 3,
                                                                       48438, 22302, 48468, 6228,
                                                                       6246, 23184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49968, 0, 3,
                                                                       48468, 22320, 48498, 6246,
                                                                       6264, 23220, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 50028, 0, 3,
                                                                       48498, 22338, 48528, 6264,
                                                                       6282, 23256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 50088, 0, 3,
                                                                       48528, 22356, 48558, 6282,
                                                                       6300, 23292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50148, 0, 3,
                                                                       48588, 22392, 48648, 6336,
                                                                       6366, 23328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50248, 0, 3,
                                                                       48648, 22428, 48708, 6366,
                                                                       6396, 23388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50348, 0, 3,
                                                                       48708, 22464, 48768, 6396,
                                                                       6426, 23448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50448, 0, 3,
                                                                       48768, 22500, 48828, 6426,
                                                                       6456, 23508, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50548, 0, 3,
                                                                       48828, 22536, 48888, 6456,
                                                                       6486, 23568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50648, 0, 3,
                                                                       48888, 22572, 48948, 6486,
                                                                       6516, 23628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50748, 0, 3,
                                                                       48948, 22608, 49008, 6516,
                                                                       6546, 23688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50848, 0, 3,
                                                                       49008, 22644, 49068, 6546,
                                                                       6576, 23748, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50948, 0, 3,
                                                                       49068, 22680, 49128, 6576,
                                                                       6606, 23808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51048, 0, 3,
                                                                       49128, 22716, 49188, 6606,
                                                                       6636, 23868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51148, 0, 3,
                                                                       49188, 22752, 49248, 6636,
                                                                       6666, 23928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51248, 0, 3,
                                                                       49248, 22788, 49308, 6666,
                                                                       6696, 23988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51348, 0, 3,
                                                                       49368, 22860, 49428, 6756,
                                                                       6786, 24048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51448, 0, 3,
                                                                       49428, 22896, 49488, 6786,
                                                                       6816, 24108, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51548, 0, 3,
                                                                       49488, 22932, 49548, 6816,
                                                                       6846, 24168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51648, 0, 3,
                                                                       49548, 22968, 49608, 6846,
                                                                       6876, 24228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51748, 0, 3,
                                                                       49608, 23004, 49668, 6876,
                                                                       6906, 24288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51848, 0, 3,
                                                                       49668, 23040, 49728, 6906,
                                                                       6936, 24348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51948, 0, 3,
                                                                       49728, 23076, 49788, 6936,
                                                                       6966, 24408, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 52048, 0, 3,
                                                                       49788, 23112, 49848, 6966,
                                                                       6996, 24468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 52148, 0, 3,
                                                                       49848, 23148, 49908, 6996,
                                                                       7026, 24528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 52248, 0, 3,
                                                                       49908, 23184, 49968, 7026,
                                                                       7056, 24588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 52348, 0, 3,
                                                                       49968, 23220, 50028, 7056,
                                                                       7086, 24648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 52448, 0, 3,
                                                                       50028, 23256, 50088, 7086,
                                                                       7116, 24708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 52548, 0, 3,
                                                                       50148, 23328, 50248, 7176,
                                                                       7221, 24768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 52698, 0, 3,
                                                                       50248, 23388, 50348, 7221,
                                                                       7266, 24858, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 52848, 0, 3,
                                                                       50348, 23448, 50448, 7266,
                                                                       7311, 24948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 52998, 0, 3,
                                                                       50448, 23508, 50548, 7311,
                                                                       7356, 25038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53148, 0, 3,
                                                                       50548, 23568, 50648, 7356,
                                                                       7401, 25128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53298, 0, 3,
                                                                       50648, 23628, 50748, 7401,
                                                                       7446, 25218, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53448, 0, 3,
                                                                       50748, 23688, 50848, 7446,
                                                                       7491, 25308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53598, 0, 3,
                                                                       50848, 23748, 50948, 7491,
                                                                       7536, 25398, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53748, 0, 3,
                                                                       50948, 23808, 51048, 7536,
                                                                       7581, 25488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53898, 0, 3,
                                                                       51048, 23868, 51148, 7581,
                                                                       7626, 25578, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54048, 0, 3,
                                                                       51148, 23928, 51248, 7626,
                                                                       7671, 25668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54198, 0, 3,
                                                                       51348, 24048, 51448, 7761,
                                                                       7806, 25758, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54348, 0, 3,
                                                                       51448, 24108, 51548, 7806,
                                                                       7851, 25848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54498, 0, 3,
                                                                       51548, 24168, 51648, 7851,
                                                                       7896, 25938, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54648, 0, 3,
                                                                       51648, 24228, 51748, 7896,
                                                                       7941, 26028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54798, 0, 3,
                                                                       51748, 24288, 51848, 7941,
                                                                       7986, 26118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54948, 0, 3,
                                                                       51848, 24348, 51948, 7986,
                                                                       8031, 26208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 55098, 0, 3,
                                                                       51948, 24408, 52048, 8031,
                                                                       8076, 26298, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 55248, 0, 3,
                                                                       52048, 24468, 52148, 8076,
                                                                       8121, 26388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 55398, 0, 3,
                                                                       52148, 24528, 52248, 8121,
                                                                       8166, 26478, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 55548, 0, 3,
                                                                       52248, 24588, 52348, 8166,
                                                                       8211, 26568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 55698, 0, 3,
                                                                       52348, 24648, 52448, 8211,
                                                                       8256, 26658, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 55848, 0, 3,
                                                                       52548, 24768, 52698, 8346,
                                                                       8409, 26748, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 56058, 0, 3,
                                                                       52698, 24858, 52848, 8409,
                                                                       8472, 26874, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 56268, 0, 3,
                                                                       52848, 24948, 52998, 8472,
                                                                       8535, 27000, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 56478, 0, 3,
                                                                       52998, 25038, 53148, 8535,
                                                                       8598, 27126, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 56688, 0, 3,
                                                                       53148, 25128, 53298, 8598,
                                                                       8661, 27252, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 56898, 0, 3,
                                                                       53298, 25218, 53448, 8661,
                                                                       8724, 27378, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 57108, 0, 3,
                                                                       53448, 25308, 53598, 8724,
                                                                       8787, 27504, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 57318, 0, 3,
                                                                       53598, 25398, 53748, 8787,
                                                                       8850, 27630, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 57528, 0, 3,
                                                                       53748, 25488, 53898, 8850,
                                                                       8913, 27756, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 57738, 0, 3,
                                                                       53898, 25578, 54048, 8913,
                                                                       8976, 27882, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 57948, 0, 3,
                                                                       54198, 25758, 54348, 9102,
                                                                       9165, 28008, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 58158, 0, 3,
                                                                       54348, 25848, 54498, 9165,
                                                                       9228, 28134, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 58368, 0, 3,
                                                                       54498, 25938, 54648, 9228,
                                                                       9291, 28260, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 58578, 0, 3,
                                                                       54648, 26028, 54798, 9291,
                                                                       9354, 28386, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 58788, 0, 3,
                                                                       54798, 26118, 54948, 9354,
                                                                       9417, 28512, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 58998, 0, 3,
                                                                       54948, 26208, 55098, 9417,
                                                                       9480, 28638, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 59208, 0, 3,
                                                                       55098, 26298, 55248, 9480,
                                                                       9543, 28764, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 59418, 0, 3,
                                                                       55248, 26388, 55398, 9543,
                                                                       9606, 28890, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 59628, 0, 3,
                                                                       55398, 26478, 55548, 9606,
                                                                       9669, 29016, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 59838, 0, 3,
                                                                       55548, 26568, 55698, 9669,
                                                                       9732, 29142, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 60048, 0, 3,
                                                                       55848, 26748, 56058, 9858,
                                                                       9942, 29268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 60328, 0, 3,
                                                                       56058, 26874, 56268, 9942,
                                                                       10026, 29436, ncols,
                                                                       gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 60608, 0, 3,
                                                                       56268, 27000, 56478,
                                                                       10026, 10110, 29604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 60888, 0, 3,
                                                                       56478, 27126, 56688,
                                                                       10110, 10194, 29772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 61168, 0, 3,
                                                                       56688, 27252, 56898,
                                                                       10194, 10278, 29940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 61448, 0, 3,
                                                                       56898, 27378, 57108,
                                                                       10278, 10362, 30108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 61728, 0, 3,
                                                                       57108, 27504, 57318,
                                                                       10362, 10446, 30276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 62008, 0, 3,
                                                                       57318, 27630, 57528,
                                                                       10446, 10530, 30444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 62288, 0, 3,
                                                                       57528, 27756, 57738,
                                                                       10530, 10614, 30612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 62568, 0, 3,
                                                                       57948, 28008, 58158,
                                                                       10782, 10866, 30780,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 62848, 0, 3,
                                                                       58158, 28134, 58368,
                                                                       10866, 10950, 30948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 63128, 0, 3,
                                                                       58368, 28260, 58578,
                                                                       10950, 11034, 31116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 63408, 0, 3,
                                                                       58578, 28386, 58788,
                                                                       11034, 11118, 31284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 63688, 0, 3,
                                                                       58788, 28512, 58998,
                                                                       11118, 11202, 31452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 63968, 0, 3,
                                                                       58998, 28638, 59208,
                                                                       11202, 11286, 31620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 64248, 0, 3,
                                                                       59208, 28764, 59418,
                                                                       11286, 11370, 31788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 64528, 0, 3,
                                                                       59418, 28890, 59628,
                                                                       11370, 11454, 31956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 64808, 0, 3,
                                                                       59628, 29016, 59838,
                                                                       11454, 11538, 32124,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 65088, 0, 3,
                                                                       60048, 29268, 60328,
                                                                       11706, 11814, 32292,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 65448, 0, 3,
                                                                       60328, 29436, 60608,
                                                                       11814, 11922, 32508,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 65808, 0, 3,
                                                                       60608, 29604, 60888,
                                                                       11922, 12030, 32724,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 66168, 0, 3,
                                                                       60888, 29772, 61168,
                                                                       12030, 12138, 32940,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 66528, 0, 3,
                                                                       61168, 29940, 61448,
                                                                       12138, 12246, 33156,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 66888, 0, 3,
                                                                       61448, 30108, 61728,
                                                                       12246, 12354, 33372,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 67248, 0, 3,
                                                                       61728, 30276, 62008,
                                                                       12354, 12462, 33588,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 67608, 0, 3,
                                                                       62008, 30444, 62288,
                                                                       12462, 12570, 33804,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 67968, 0, 3,
                                                                       62568, 30780, 62848,
                                                                       12786, 12894, 34020,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 68328, 0, 3,
                                                                       62848, 30948, 63128,
                                                                       12894, 13002, 34236,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 68688, 0, 3,
                                                                       63128, 31116, 63408,
                                                                       13002, 13110, 34452,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 69048, 0, 3,
                                                                       63408, 31284, 63688,
                                                                       13110, 13218, 34668,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 69408, 0, 3,
                                                                       63688, 31452, 63968,
                                                                       13218, 13326, 34884,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 69768, 0, 3,
                                                                       63968, 31620, 64248,
                                                                       13326, 13434, 35100,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 70128, 0, 3,
                                                                       64248, 31788, 64528,
                                                                       13434, 13542, 35316,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 70488, 0, 3,
                                                                       64528, 31956, 64808,
                                                                       13542, 13650, 35532,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 70848, 0, 3,
                                                                       65088, 32292, 65448,
                                                                       13866, 14001, 35748,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 71298, 0, 3,
                                                                       65448, 32508, 65808,
                                                                       14001, 14136, 36018,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 71748, 0, 3,
                                                                       65808, 32724, 66168,
                                                                       14136, 14271, 36288,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 72198, 0, 3,
                                                                       66168, 32940, 66528,
                                                                       14271, 14406, 36558,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 72648, 0, 3,
                                                                       66528, 33156, 66888,
                                                                       14406, 14541, 36828,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 73098, 0, 3,
                                                                       66888, 33372, 67248,
                                                                       14541, 14676, 37098,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 73548, 0, 3,
                                                                       67248, 33588, 67608,
                                                                       14676, 14811, 37368,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 73998, 0, 3,
                                                                       67968, 34020, 68328,
                                                                       15081, 15216, 37638,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 74448, 0, 3,
                                                                       68328, 34236, 68688,
                                                                       15216, 15351, 37908,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 74898, 0, 3,
                                                                       68688, 34452, 69048,
                                                                       15351, 15486, 38178,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 75348, 0, 3,
                                                                       69048, 34668, 69408,
                                                                       15486, 15621, 38448,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 75798, 0, 3,
                                                                       69408, 34884, 69768,
                                                                       15621, 15756, 38718,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 76248, 0, 3,
                                                                       69768, 35100, 70128,
                                                                       15756, 15891, 38988,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 76698, 0, 3,
                                                                       70128, 35316, 70488,
                                                                       15891, 16026, 39258,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 77148, 0, 3,
                                                                       70848, 35748, 71298,
                                                                       16296, 16461, 39528,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 77698, 0, 3,
                                                                       71298, 36018, 71748,
                                                                       16461, 16626, 39858,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 78248, 0, 3,
                                                                       71748, 36288, 72198,
                                                                       16626, 16791, 40188,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 78798, 0, 3,
                                                                       72198, 36558, 72648,
                                                                       16791, 16956, 40518,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 79348, 0, 3,
                                                                       72648, 36828, 73098,
                                                                       16956, 17121, 40848,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 79898, 0, 3,
                                                                       73098, 37098, 73548,
                                                                       17121, 17286, 41178,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 80448, 0, 3,
                                                                       73998, 37638, 74448,
                                                                       17616, 17781, 41508,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 80998, 0, 3,
                                                                       74448, 37908, 74898,
                                                                       17781, 17946, 41838,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 81548, 0, 3,
                                                                       74898, 38178, 75348,
                                                                       17946, 18111, 42168,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 82098, 0, 3,
                                                                       75348, 38448, 75798,
                                                                       18111, 18276, 42498,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 82648, 0, 3,
                                                                       75798, 38718, 76248,
                                                                       18276, 18441, 42828,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 83198, 0, 3,
                                                                       76248, 38988, 76698,
                                                                       18441, 18606, 43158,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 83748, 0, 3,
                                                                       77148, 39528, 77698,
                                                                       18936, 19134, 43488,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 84408, 0, 3,
                                                                       77698, 39858, 78248,
                                                                       19134, 19332, 43884,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 85068, 0, 3,
                                                                       78248, 40188, 78798,
                                                                       19332, 19530, 44280,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 85728, 0, 3,
                                                                       78798, 40518, 79348,
                                                                       19530, 19728, 44676,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 86388, 0, 3,
                                                                       79348, 40848, 79898,
                                                                       19728, 19926, 45072,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 87048, 0, 3,
                                                                       80448, 41508, 80998,
                                                                       20322, 20520, 45468,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 87708, 0, 3,
                                                                       80998, 41838, 81548,
                                                                       20520, 20718, 45864,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 88368, 0, 3,
                                                                       81548, 42168, 82098,
                                                                       20718, 20916, 46260,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 89028, 0, 3,
                                                                       82098, 42498, 82648,
                                                                       20916, 21114, 46656,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 89688, 0, 3,
                                                                       82648, 42828, 83198,
                                                                       21114, 21312, 47052,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90348, 3, 21708,
                                                                       21714, 47468, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90363, 3, 21714,
                                                                       21720, 47478, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90378, 3, 21720,
                                                                       21726, 47488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90393, 3, 21726,
                                                                       21732, 47498, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90408, 3, 21732,
                                                                       21738, 47508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90423, 3, 21738,
                                                                       21744, 47518, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90438, 3, 21744,
                                                                       21750, 47528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90453, 3, 21750,
                                                                       21756, 47538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90468, 3, 21756,
                                                                       21762, 47548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90483, 3, 21762,
                                                                       21768, 47558, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90498, 3, 21768,
                                                                       21774, 47568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90513, 3, 21774,
                                                                       21780, 47578, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90528, 3, 21780,
                                                                       21786, 47588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90543, 3, 21798,
                                                                       21804, 47618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90558, 3, 21804,
                                                                       21810, 47628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90573, 3, 21810,
                                                                       21816, 47638, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90588, 3, 21816,
                                                                       21822, 47648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90603, 3, 21822,
                                                                       21828, 47658, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90618, 3, 21828,
                                                                       21834, 47668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90633, 3, 21834,
                                                                       21840, 47678, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90648, 3, 21840,
                                                                       21846, 47688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90663, 3, 21846,
                                                                       21852, 47698, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90678, 3, 21852,
                                                                       21858, 47708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90693, 3, 21858,
                                                                       21864, 47718, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90708, 3, 21864,
                                                                       21870, 47728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 90723, 3, 21870,
                                                                       21876, 47738, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90738, 0, 3,
                                                                       90348, 47468, 90363,
                                                                       21888, 21906, 47808,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90783, 0, 3,
                                                                       90363, 47478, 90378,
                                                                       21906, 21924, 47838,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90828, 0, 3,
                                                                       90378, 47488, 90393,
                                                                       21924, 21942, 47868,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90873, 0, 3,
                                                                       90393, 47498, 90408,
                                                                       21942, 21960, 47898,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90918, 0, 3,
                                                                       90408, 47508, 90423,
                                                                       21960, 21978, 47928,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90963, 0, 3,
                                                                       90423, 47518, 90438,
                                                                       21978, 21996, 47958,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91008, 0, 3,
                                                                       90438, 47528, 90453,
                                                                       21996, 22014, 47988,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91053, 0, 3,
                                                                       90453, 47538, 90468,
                                                                       22014, 22032, 48018,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91098, 0, 3,
                                                                       90468, 47548, 90483,
                                                                       22032, 22050, 48048,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91143, 0, 3,
                                                                       90483, 47558, 90498,
                                                                       22050, 22068, 48078,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91188, 0, 3,
                                                                       90498, 47568, 90513,
                                                                       22068, 22086, 48108,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91233, 0, 3,
                                                                       90513, 47578, 90528,
                                                                       22086, 22104, 48138,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91278, 0, 3,
                                                                       90543, 47618, 90558,
                                                                       22140, 22158, 48228,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91323, 0, 3,
                                                                       90558, 47628, 90573,
                                                                       22158, 22176, 48258,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91368, 0, 3,
                                                                       90573, 47638, 90588,
                                                                       22176, 22194, 48288,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91413, 0, 3,
                                                                       90588, 47648, 90603,
                                                                       22194, 22212, 48318,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91458, 0, 3,
                                                                       90603, 47658, 90618,
                                                                       22212, 22230, 48348,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91503, 0, 3,
                                                                       90618, 47668, 90633,
                                                                       22230, 22248, 48378,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91548, 0, 3,
                                                                       90633, 47678, 90648,
                                                                       22248, 22266, 48408,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91593, 0, 3,
                                                                       90648, 47688, 90663,
                                                                       22266, 22284, 48438,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91638, 0, 3,
                                                                       90663, 47698, 90678,
                                                                       22284, 22302, 48468,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91683, 0, 3,
                                                                       90678, 47708, 90693,
                                                                       22302, 22320, 48498,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91728, 0, 3,
                                                                       90693, 47718, 90708,
                                                                       22320, 22338, 48528,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 91773, 0, 3,
                                                                       90708, 47728, 90723,
                                                                       22338, 22356, 48558,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91818, 0, 3,
                                                                       90738, 47808, 90783,
                                                                       22392, 22428, 48708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91908, 0, 3,
                                                                       90783, 47838, 90828,
                                                                       22428, 22464, 48768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91998, 0, 3,
                                                                       90828, 47868, 90873,
                                                                       22464, 22500, 48828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92088, 0, 3,
                                                                       90873, 47898, 90918,
                                                                       22500, 22536, 48888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92178, 0, 3,
                                                                       90918, 47928, 90963,
                                                                       22536, 22572, 48948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92268, 0, 3,
                                                                       90963, 47958, 91008,
                                                                       22572, 22608, 49008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92358, 0, 3,
                                                                       91008, 47988, 91053,
                                                                       22608, 22644, 49068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92448, 0, 3,
                                                                       91053, 48018, 91098,
                                                                       22644, 22680, 49128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92538, 0, 3,
                                                                       91098, 48048, 91143,
                                                                       22680, 22716, 49188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92628, 0, 3,
                                                                       91143, 48078, 91188,
                                                                       22716, 22752, 49248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92718, 0, 3,
                                                                       91188, 48108, 91233,
                                                                       22752, 22788, 49308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92808, 0, 3,
                                                                       91278, 48228, 91323,
                                                                       22860, 22896, 49488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92898, 0, 3,
                                                                       91323, 48258, 91368,
                                                                       22896, 22932, 49548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92988, 0, 3,
                                                                       91368, 48288, 91413,
                                                                       22932, 22968, 49608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 93078, 0, 3,
                                                                       91413, 48318, 91458,
                                                                       22968, 23004, 49668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 93168, 0, 3,
                                                                       91458, 48348, 91503,
                                                                       23004, 23040, 49728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 93258, 0, 3,
                                                                       91503, 48378, 91548,
                                                                       23040, 23076, 49788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 93348, 0, 3,
                                                                       91548, 48408, 91593,
                                                                       23076, 23112, 49848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 93438, 0, 3,
                                                                       91593, 48438, 91638,
                                                                       23112, 23148, 49908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 93528, 0, 3,
                                                                       91638, 48468, 91683,
                                                                       23148, 23184, 49968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 93618, 0, 3,
                                                                       91683, 48498, 91728,
                                                                       23184, 23220, 50028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 93708, 0, 3,
                                                                       91728, 48528, 91773,
                                                                       23220, 23256, 50088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 93798, 0, 3,
                                                                       91818, 48708, 91908,
                                                                       23328, 23388, 50348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 93948, 0, 3,
                                                                       91908, 48768, 91998,
                                                                       23388, 23448, 50448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94098, 0, 3,
                                                                       91998, 48828, 92088,
                                                                       23448, 23508, 50548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94248, 0, 3,
                                                                       92088, 48888, 92178,
                                                                       23508, 23568, 50648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94398, 0, 3,
                                                                       92178, 48948, 92268,
                                                                       23568, 23628, 50748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94548, 0, 3,
                                                                       92268, 49008, 92358,
                                                                       23628, 23688, 50848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94698, 0, 3,
                                                                       92358, 49068, 92448,
                                                                       23688, 23748, 50948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94848, 0, 3,
                                                                       92448, 49128, 92538,
                                                                       23748, 23808, 51048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94998, 0, 3,
                                                                       92538, 49188, 92628,
                                                                       23808, 23868, 51148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 95148, 0, 3,
                                                                       92628, 49248, 92718,
                                                                       23868, 23928, 51248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 95298, 0, 3,
                                                                       92808, 49488, 92898,
                                                                       24048, 24108, 51548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 95448, 0, 3,
                                                                       92898, 49548, 92988,
                                                                       24108, 24168, 51648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 95598, 0, 3,
                                                                       92988, 49608, 93078,
                                                                       24168, 24228, 51748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 95748, 0, 3,
                                                                       93078, 49668, 93168,
                                                                       24228, 24288, 51848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 95898, 0, 3,
                                                                       93168, 49728, 93258,
                                                                       24288, 24348, 51948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 96048, 0, 3,
                                                                       93258, 49788, 93348,
                                                                       24348, 24408, 52048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 96198, 0, 3,
                                                                       93348, 49848, 93438,
                                                                       24408, 24468, 52148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 96348, 0, 3,
                                                                       93438, 49908, 93528,
                                                                       24468, 24528, 52248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 96498, 0, 3,
                                                                       93528, 49968, 93618,
                                                                       24528, 24588, 52348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 96648, 0, 3,
                                                                       93618, 50028, 93708,
                                                                       24588, 24648, 52448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 96798, 0, 3,
                                                                       93798, 50348, 93948,
                                                                       24768, 24858, 52848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 97023, 0, 3,
                                                                       93948, 50448, 94098,
                                                                       24858, 24948, 52998,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 97248, 0, 3,
                                                                       94098, 50548, 94248,
                                                                       24948, 25038, 53148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 97473, 0, 3,
                                                                       94248, 50648, 94398,
                                                                       25038, 25128, 53298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 97698, 0, 3,
                                                                       94398, 50748, 94548,
                                                                       25128, 25218, 53448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 97923, 0, 3,
                                                                       94548, 50848, 94698,
                                                                       25218, 25308, 53598,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 98148, 0, 3,
                                                                       94698, 50948, 94848,
                                                                       25308, 25398, 53748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 98373, 0, 3,
                                                                       94848, 51048, 94998,
                                                                       25398, 25488, 53898,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 98598, 0, 3,
                                                                       94998, 51148, 95148,
                                                                       25488, 25578, 54048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 98823, 0, 3,
                                                                       95298, 51548, 95448,
                                                                       25758, 25848, 54498,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 99048, 0, 3,
                                                                       95448, 51648, 95598,
                                                                       25848, 25938, 54648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 99273, 0, 3,
                                                                       95598, 51748, 95748,
                                                                       25938, 26028, 54798,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 99498, 0, 3,
                                                                       95748, 51848, 95898,
                                                                       26028, 26118, 54948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 99723, 0, 3,
                                                                       95898, 51948, 96048,
                                                                       26118, 26208, 55098,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 99948, 0, 3,
                                                                       96048, 52048, 96198,
                                                                       26208, 26298, 55248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 100173, 0, 3,
                                                                       96198, 52148, 96348,
                                                                       26298, 26388, 55398,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 100398, 0, 3,
                                                                       96348, 52248, 96498,
                                                                       26388, 26478, 55548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 100623, 0, 3,
                                                                       96498, 52348, 96648,
                                                                       26478, 26568, 55698,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 100848, 0, 3,
                                                                       96798, 52848, 97023,
                                                                       26748, 26874, 56268,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 101163, 0, 3,
                                                                       97023, 52998, 97248,
                                                                       26874, 27000, 56478,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 101478, 0, 3,
                                                                       97248, 53148, 97473,
                                                                       27000, 27126, 56688,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 101793, 0, 3,
                                                                       97473, 53298, 97698,
                                                                       27126, 27252, 56898,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 102108, 0, 3,
                                                                       97698, 53448, 97923,
                                                                       27252, 27378, 57108,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 102423, 0, 3,
                                                                       97923, 53598, 98148,
                                                                       27378, 27504, 57318,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 102738, 0, 3,
                                                                       98148, 53748, 98373,
                                                                       27504, 27630, 57528,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 103053, 0, 3,
                                                                       98373, 53898, 98598,
                                                                       27630, 27756, 57738,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 103368, 0, 3,
                                                                       98823, 54498, 99048,
                                                                       28008, 28134, 58368,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 103683, 0, 3,
                                                                       99048, 54648, 99273,
                                                                       28134, 28260, 58578,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 103998, 0, 3,
                                                                       99273, 54798, 99498,
                                                                       28260, 28386, 58788,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 104313, 0, 3,
                                                                       99498, 54948, 99723,
                                                                       28386, 28512, 58998,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 104628, 0, 3,
                                                                       99723, 55098, 99948,
                                                                       28512, 28638, 59208,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 104943, 0, 3,
                                                                       99948, 55248, 100173,
                                                                       28638, 28764, 59418,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 105258, 0, 3,
                                                                       100173, 55398, 100398,
                                                                       28764, 28890, 59628,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 105573, 0, 3,
                                                                       100398, 55548, 100623,
                                                                       28890, 29016, 59838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 105888, 0, 3,
                                                                       100848, 56268, 101163,
                                                                       29268, 29436, 60608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 106308, 0, 3,
                                                                       101163, 56478, 101478,
                                                                       29436, 29604, 60888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 106728, 0, 3,
                                                                       101478, 56688, 101793,
                                                                       29604, 29772, 61168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 107148, 0, 3,
                                                                       101793, 56898, 102108,
                                                                       29772, 29940, 61448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 107568, 0, 3,
                                                                       102108, 57108, 102423,
                                                                       29940, 30108, 61728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 107988, 0, 3,
                                                                       102423, 57318, 102738,
                                                                       30108, 30276, 62008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 108408, 0, 3,
                                                                       102738, 57528, 103053,
                                                                       30276, 30444, 62288,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 108828, 0, 3,
                                                                       103368, 58368, 103683,
                                                                       30780, 30948, 63128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 109248, 0, 3,
                                                                       103683, 58578, 103998,
                                                                       30948, 31116, 63408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 109668, 0, 3,
                                                                       103998, 58788, 104313,
                                                                       31116, 31284, 63688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 110088, 0, 3,
                                                                       104313, 58998, 104628,
                                                                       31284, 31452, 63968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 110508, 0, 3,
                                                                       104628, 59208, 104943,
                                                                       31452, 31620, 64248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 110928, 0, 3,
                                                                       104943, 59418, 105258,
                                                                       31620, 31788, 64528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 111348, 0, 3,
                                                                       105258, 59628, 105573,
                                                                       31788, 31956, 64808,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 111768, 0, 3,
                                                                       105888, 60608, 106308,
                                                                       32292, 32508, 65808,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 112308, 0, 3,
                                                                       106308, 60888, 106728,
                                                                       32508, 32724, 66168,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 112848, 0, 3,
                                                                       106728, 61168, 107148,
                                                                       32724, 32940, 66528,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 113388, 0, 3,
                                                                       107148, 61448, 107568,
                                                                       32940, 33156, 66888,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 113928, 0, 3,
                                                                       107568, 61728, 107988,
                                                                       33156, 33372, 67248,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 114468, 0, 3,
                                                                       107988, 62008, 108408,
                                                                       33372, 33588, 67608,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 115008, 0, 3,
                                                                       108828, 63128, 109248,
                                                                       34020, 34236, 68688,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 115548, 0, 3,
                                                                       109248, 63408, 109668,
                                                                       34236, 34452, 69048,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 116088, 0, 3,
                                                                       109668, 63688, 110088,
                                                                       34452, 34668, 69408,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 116628, 0, 3,
                                                                       110088, 63968, 110508,
                                                                       34668, 34884, 69768,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 117168, 0, 3,
                                                                       110508, 64248, 110928,
                                                                       34884, 35100, 70128,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 117708, 0, 3,
                                                                       110928, 64528, 111348,
                                                                       35100, 35316, 70488,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 118248, 0, 3,
                                                                       111768, 65808, 112308,
                                                                       35748, 36018, 71748,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 118923, 0, 3,
                                                                       112308, 66168, 112848,
                                                                       36018, 36288, 72198,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 119598, 0, 3,
                                                                       112848, 66528, 113388,
                                                                       36288, 36558, 72648,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 120273, 0, 3,
                                                                       113388, 66888, 113928,
                                                                       36558, 36828, 73098,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 120948, 0, 3,
                                                                       113928, 67248, 114468,
                                                                       36828, 37098, 73548,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 121623, 0, 3,
                                                                       115008, 68688, 115548,
                                                                       37638, 37908, 74898,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 122298, 0, 3,
                                                                       115548, 69048, 116088,
                                                                       37908, 38178, 75348,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 122973, 0, 3,
                                                                       116088, 69408, 116628,
                                                                       38178, 38448, 75798,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 123648, 0, 3,
                                                                       116628, 69768, 117168,
                                                                       38448, 38718, 76248,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 124323, 0, 3,
                                                                       117168, 70128, 117708,
                                                                       38718, 38988, 76698,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 124998, 0, 3,
                                                                       118248, 71748, 118923,
                                                                       39528, 39858, 78248,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 125823, 0, 3,
                                                                       118923, 72198, 119598,
                                                                       39858, 40188, 78798,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 126648, 0, 3,
                                                                       119598, 72648, 120273,
                                                                       40188, 40518, 79348,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 127473, 0, 3,
                                                                       120273, 73098, 120948,
                                                                       40518, 40848, 79898,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 128298, 0, 3,
                                                                       121623, 74898, 122298,
                                                                       41508, 41838, 81548,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 129123, 0, 3,
                                                                       122298, 75348, 122973,
                                                                       41838, 42168, 82098,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 129948, 0, 3,
                                                                       122973, 75798, 123648,
                                                                       42168, 42498, 82648,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 130773, 0, 3,
                                                                       123648, 76248, 124323,
                                                                       42498, 42828, 83198,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 131598, 0, 3,
                                                                       124998, 78248, 125823,
                                                                       43488, 43884, 85068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 132588, 0, 3,
                                                                       125823, 78798, 126648,
                                                                       43884, 44280, 85728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 133578, 0, 3,
                                                                       126648, 79348, 127473,
                                                                       44280, 44676, 86388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 134568, 0, 3,
                                                                       128298, 81548, 129123,
                                                                       45468, 45864, 88368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 135558, 0, 3,
                                                                       129123, 82098, 129948,
                                                                       45864, 46260, 89028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 136548, 0, 3,
                                                                       129948, 82648, 130773,
                                                                       46260, 46656, 89688,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137538, 3, 47448,
                                                                       47458, 90348, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137559, 3, 47458,
                                                                       47468, 90363, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137580, 3, 47468,
                                                                       47478, 90378, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137601, 3, 47478,
                                                                       47488, 90393, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137622, 3, 47488,
                                                                       47498, 90408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137643, 3, 47498,
                                                                       47508, 90423, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137664, 3, 47508,
                                                                       47518, 90438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137685, 3, 47518,
                                                                       47528, 90453, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137706, 3, 47528,
                                                                       47538, 90468, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137727, 3, 47538,
                                                                       47548, 90483, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137748, 3, 47548,
                                                                       47558, 90498, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137769, 3, 47558,
                                                                       47568, 90513, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137790, 3, 47568,
                                                                       47578, 90528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137811, 3, 47598,
                                                                       47608, 90543, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137832, 3, 47608,
                                                                       47618, 90558, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137853, 3, 47618,
                                                                       47628, 90573, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137874, 3, 47628,
                                                                       47638, 90588, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137895, 3, 47638,
                                                                       47648, 90603, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137916, 3, 47648,
                                                                       47658, 90618, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137937, 3, 47658,
                                                                       47668, 90633, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137958, 3, 47668,
                                                                       47678, 90648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 137979, 3, 47678,
                                                                       47688, 90663, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 138000, 3, 47688,
                                                                       47698, 90678, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 138021, 3, 47698,
                                                                       47708, 90693, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 138042, 3, 47708,
                                                                       47718, 90708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 138063, 3, 47718,
                                                                       47728, 90723, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138084, 0, 3,
                                                                       137538, 90348, 137559,
                                                                       47748, 47778, 90738,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138147, 0, 3,
                                                                       137559, 90363, 137580,
                                                                       47778, 47808, 90783,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138210, 0, 3,
                                                                       137580, 90378, 137601,
                                                                       47808, 47838, 90828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138273, 0, 3,
                                                                       137601, 90393, 137622,
                                                                       47838, 47868, 90873,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138336, 0, 3,
                                                                       137622, 90408, 137643,
                                                                       47868, 47898, 90918,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138399, 0, 3,
                                                                       137643, 90423, 137664,
                                                                       47898, 47928, 90963,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138462, 0, 3,
                                                                       137664, 90438, 137685,
                                                                       47928, 47958, 91008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138525, 0, 3,
                                                                       137685, 90453, 137706,
                                                                       47958, 47988, 91053,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138588, 0, 3,
                                                                       137706, 90468, 137727,
                                                                       47988, 48018, 91098,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138651, 0, 3,
                                                                       137727, 90483, 137748,
                                                                       48018, 48048, 91143,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138714, 0, 3,
                                                                       137748, 90498, 137769,
                                                                       48048, 48078, 91188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138777, 0, 3,
                                                                       137769, 90513, 137790,
                                                                       48078, 48108, 91233,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138840, 0, 3,
                                                                       137811, 90543, 137832,
                                                                       48168, 48198, 91278,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138903, 0, 3,
                                                                       137832, 90558, 137853,
                                                                       48198, 48228, 91323,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 138966, 0, 3,
                                                                       137853, 90573, 137874,
                                                                       48228, 48258, 91368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 139029, 0, 3,
                                                                       137874, 90588, 137895,
                                                                       48258, 48288, 91413,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 139092, 0, 3,
                                                                       137895, 90603, 137916,
                                                                       48288, 48318, 91458,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 139155, 0, 3,
                                                                       137916, 90618, 137937,
                                                                       48318, 48348, 91503,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 139218, 0, 3,
                                                                       137937, 90633, 137958,
                                                                       48348, 48378, 91548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 139281, 0, 3,
                                                                       137958, 90648, 137979,
                                                                       48378, 48408, 91593,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 139344, 0, 3,
                                                                       137979, 90663, 138000,
                                                                       48408, 48438, 91638,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 139407, 0, 3,
                                                                       138000, 90678, 138021,
                                                                       48438, 48468, 91683,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 139470, 0, 3,
                                                                       138021, 90693, 138042,
                                                                       48468, 48498, 91728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 139533, 0, 3,
                                                                       138042, 90708, 138063,
                                                                       48498, 48528, 91773,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 139596, 0, 3,
                                                                       138084, 90738, 138147,
                                                                       48588, 48648, 91818,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 139722, 0, 3,
                                                                       138147, 90783, 138210,
                                                                       48648, 48708, 91908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 139848, 0, 3,
                                                                       138210, 90828, 138273,
                                                                       48708, 48768, 91998,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 139974, 0, 3,
                                                                       138273, 90873, 138336,
                                                                       48768, 48828, 92088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 140100, 0, 3,
                                                                       138336, 90918, 138399,
                                                                       48828, 48888, 92178,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 140226, 0, 3,
                                                                       138399, 90963, 138462,
                                                                       48888, 48948, 92268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 140352, 0, 3,
                                                                       138462, 91008, 138525,
                                                                       48948, 49008, 92358,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 140478, 0, 3,
                                                                       138525, 91053, 138588,
                                                                       49008, 49068, 92448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 140604, 0, 3,
                                                                       138588, 91098, 138651,
                                                                       49068, 49128, 92538,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 140730, 0, 3,
                                                                       138651, 91143, 138714,
                                                                       49128, 49188, 92628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 140856, 0, 3,
                                                                       138714, 91188, 138777,
                                                                       49188, 49248, 92718,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 140982, 0, 3,
                                                                       138840, 91278, 138903,
                                                                       49368, 49428, 92808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 141108, 0, 3,
                                                                       138903, 91323, 138966,
                                                                       49428, 49488, 92898,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 141234, 0, 3,
                                                                       138966, 91368, 139029,
                                                                       49488, 49548, 92988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 141360, 0, 3,
                                                                       139029, 91413, 139092,
                                                                       49548, 49608, 93078,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 141486, 0, 3,
                                                                       139092, 91458, 139155,
                                                                       49608, 49668, 93168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 141612, 0, 3,
                                                                       139155, 91503, 139218,
                                                                       49668, 49728, 93258,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 141738, 0, 3,
                                                                       139218, 91548, 139281,
                                                                       49728, 49788, 93348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 141864, 0, 3,
                                                                       139281, 91593, 139344,
                                                                       49788, 49848, 93438,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 141990, 0, 3,
                                                                       139344, 91638, 139407,
                                                                       49848, 49908, 93528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 142116, 0, 3,
                                                                       139407, 91683, 139470,
                                                                       49908, 49968, 93618,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 142242, 0, 3,
                                                                       139470, 91728, 139533,
                                                                       49968, 50028, 93708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 142368, 0, 3,
                                                                       139596, 91818, 139722,
                                                                       50148, 50248, 93798,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 142578, 0, 3,
                                                                       139722, 91908, 139848,
                                                                       50248, 50348, 93948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 142788, 0, 3,
                                                                       139848, 91998, 139974,
                                                                       50348, 50448, 94098,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 142998, 0, 3,
                                                                       139974, 92088, 140100,
                                                                       50448, 50548, 94248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 143208, 0, 3,
                                                                       140100, 92178, 140226,
                                                                       50548, 50648, 94398,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 143418, 0, 3,
                                                                       140226, 92268, 140352,
                                                                       50648, 50748, 94548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 143628, 0, 3,
                                                                       140352, 92358, 140478,
                                                                       50748, 50848, 94698,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 143838, 0, 3,
                                                                       140478, 92448, 140604,
                                                                       50848, 50948, 94848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 144048, 0, 3,
                                                                       140604, 92538, 140730,
                                                                       50948, 51048, 94998,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 144258, 0, 3,
                                                                       140730, 92628, 140856,
                                                                       51048, 51148, 95148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 144468, 0, 3,
                                                                       140982, 92808, 141108,
                                                                       51348, 51448, 95298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 144678, 0, 3,
                                                                       141108, 92898, 141234,
                                                                       51448, 51548, 95448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 144888, 0, 3,
                                                                       141234, 92988, 141360,
                                                                       51548, 51648, 95598,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 145098, 0, 3,
                                                                       141360, 93078, 141486,
                                                                       51648, 51748, 95748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 145308, 0, 3,
                                                                       141486, 93168, 141612,
                                                                       51748, 51848, 95898,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 145518, 0, 3,
                                                                       141612, 93258, 141738,
                                                                       51848, 51948, 96048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 145728, 0, 3,
                                                                       141738, 93348, 141864,
                                                                       51948, 52048, 96198,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 145938, 0, 3,
                                                                       141864, 93438, 141990,
                                                                       52048, 52148, 96348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 146148, 0, 3,
                                                                       141990, 93528, 142116,
                                                                       52148, 52248, 96498,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 146358, 0, 3,
                                                                       142116, 93618, 142242,
                                                                       52248, 52348, 96648,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 146568, 0, 3,
                                                                       142368, 93798, 142578,
                                                                       52548, 52698, 96798,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 146883, 0, 3,
                                                                       142578, 93948, 142788,
                                                                       52698, 52848, 97023,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 147198, 0, 3,
                                                                       142788, 94098, 142998,
                                                                       52848, 52998, 97248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 147513, 0, 3,
                                                                       142998, 94248, 143208,
                                                                       52998, 53148, 97473,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 147828, 0, 3,
                                                                       143208, 94398, 143418,
                                                                       53148, 53298, 97698,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 148143, 0, 3,
                                                                       143418, 94548, 143628,
                                                                       53298, 53448, 97923,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 148458, 0, 3,
                                                                       143628, 94698, 143838,
                                                                       53448, 53598, 98148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 148773, 0, 3,
                                                                       143838, 94848, 144048,
                                                                       53598, 53748, 98373,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 149088, 0, 3,
                                                                       144048, 94998, 144258,
                                                                       53748, 53898, 98598,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 149403, 0, 3,
                                                                       144468, 95298, 144678,
                                                                       54198, 54348, 98823,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 149718, 0, 3,
                                                                       144678, 95448, 144888,
                                                                       54348, 54498, 99048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 150033, 0, 3,
                                                                       144888, 95598, 145098,
                                                                       54498, 54648, 99273,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 150348, 0, 3,
                                                                       145098, 95748, 145308,
                                                                       54648, 54798, 99498,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 150663, 0, 3,
                                                                       145308, 95898, 145518,
                                                                       54798, 54948, 99723,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 150978, 0, 3,
                                                                       145518, 96048, 145728,
                                                                       54948, 55098, 99948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 151293, 0, 3,
                                                                       145728, 96198, 145938,
                                                                       55098, 55248, 100173,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 151608, 0, 3,
                                                                       145938, 96348, 146148,
                                                                       55248, 55398, 100398,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 151923, 0, 3,
                                                                       146148, 96498, 146358,
                                                                       55398, 55548, 100623,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 152238, 0, 3,
                                                                       146568, 96798, 146883,
                                                                       55848, 56058, 100848,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 152679, 0, 3,
                                                                       146883, 97023, 147198,
                                                                       56058, 56268, 101163,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 153120, 0, 3,
                                                                       147198, 97248, 147513,
                                                                       56268, 56478, 101478,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 153561, 0, 3,
                                                                       147513, 97473, 147828,
                                                                       56478, 56688, 101793,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 154002, 0, 3,
                                                                       147828, 97698, 148143,
                                                                       56688, 56898, 102108,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 154443, 0, 3,
                                                                       148143, 97923, 148458,
                                                                       56898, 57108, 102423,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 154884, 0, 3,
                                                                       148458, 98148, 148773,
                                                                       57108, 57318, 102738,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 155325, 0, 3,
                                                                       148773, 98373, 149088,
                                                                       57318, 57528, 103053,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 155766, 0, 3,
                                                                       149403, 98823, 149718,
                                                                       57948, 58158, 103368,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 156207, 0, 3,
                                                                       149718, 99048, 150033,
                                                                       58158, 58368, 103683,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 156648, 0, 3,
                                                                       150033, 99273, 150348,
                                                                       58368, 58578, 103998,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 157089, 0, 3,
                                                                       150348, 99498, 150663,
                                                                       58578, 58788, 104313,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 157530, 0, 3,
                                                                       150663, 99723, 150978,
                                                                       58788, 58998, 104628,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 157971, 0, 3,
                                                                       150978, 99948, 151293,
                                                                       58998, 59208, 104943,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 158412, 0, 3,
                                                                       151293, 100173, 151608,
                                                                       59208, 59418, 105258,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 158853, 0, 3,
                                                                       151608, 100398, 151923,
                                                                       59418, 59628, 105573,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 159294, 0, 3,
                                                                       152238, 100848, 152679,
                                                                       60048, 60328, 105888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 159882, 0, 3,
                                                                       152679, 101163, 153120,
                                                                       60328, 60608, 106308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 160470, 0, 3,
                                                                       153120, 101478, 153561,
                                                                       60608, 60888, 106728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 161058, 0, 3,
                                                                       153561, 101793, 154002,
                                                                       60888, 61168, 107148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 161646, 0, 3,
                                                                       154002, 102108, 154443,
                                                                       61168, 61448, 107568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 162234, 0, 3,
                                                                       154443, 102423, 154884,
                                                                       61448, 61728, 107988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 162822, 0, 3,
                                                                       154884, 102738, 155325,
                                                                       61728, 62008, 108408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 163410, 0, 3,
                                                                       155766, 103368, 156207,
                                                                       62568, 62848, 108828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 163998, 0, 3,
                                                                       156207, 103683, 156648,
                                                                       62848, 63128, 109248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 164586, 0, 3,
                                                                       156648, 103998, 157089,
                                                                       63128, 63408, 109668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 165174, 0, 3,
                                                                       157089, 104313, 157530,
                                                                       63408, 63688, 110088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 165762, 0, 3,
                                                                       157530, 104628, 157971,
                                                                       63688, 63968, 110508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 166350, 0, 3,
                                                                       157971, 104943, 158412,
                                                                       63968, 64248, 110928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 166938, 0, 3,
                                                                       158412, 105258, 158853,
                                                                       64248, 64528, 111348,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 167526, 0, 3,
                                                                       159294, 105888, 159882,
                                                                       65088, 65448, 111768,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 168282, 0, 3,
                                                                       159882, 106308, 160470,
                                                                       65448, 65808, 112308,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 169038, 0, 3,
                                                                       160470, 106728, 161058,
                                                                       65808, 66168, 112848,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 169794, 0, 3,
                                                                       161058, 107148, 161646,
                                                                       66168, 66528, 113388,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 170550, 0, 3,
                                                                       161646, 107568, 162234,
                                                                       66528, 66888, 113928,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 171306, 0, 3,
                                                                       162234, 107988, 162822,
                                                                       66888, 67248, 114468,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 172062, 0, 3,
                                                                       163410, 108828, 163998,
                                                                       67968, 68328, 115008,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 172818, 0, 3,
                                                                       163998, 109248, 164586,
                                                                       68328, 68688, 115548,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 173574, 0, 3,
                                                                       164586, 109668, 165174,
                                                                       68688, 69048, 116088,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 174330, 0, 3,
                                                                       165174, 110088, 165762,
                                                                       69048, 69408, 116628,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 175086, 0, 3,
                                                                       165762, 110508, 166350,
                                                                       69408, 69768, 117168,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 175842, 0, 3,
                                                                       166350, 110928, 166938,
                                                                       69768, 70128, 117708,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 176598, 0, 3,
                                                                       167526, 111768, 168282,
                                                                       70848, 71298, 118248,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 177543, 0, 3,
                                                                       168282, 112308, 169038,
                                                                       71298, 71748, 118923,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 178488, 0, 3,
                                                                       169038, 112848, 169794,
                                                                       71748, 72198, 119598,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 179433, 0, 3,
                                                                       169794, 113388, 170550,
                                                                       72198, 72648, 120273,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 180378, 0, 3,
                                                                       170550, 113928, 171306,
                                                                       72648, 73098, 120948,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 181323, 0, 3,
                                                                       172062, 115008, 172818,
                                                                       73998, 74448, 121623,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 182268, 0, 3,
                                                                       172818, 115548, 173574,
                                                                       74448, 74898, 122298,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 183213, 0, 3,
                                                                       173574, 116088, 174330,
                                                                       74898, 75348, 122973,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 184158, 0, 3,
                                                                       174330, 116628, 175086,
                                                                       75348, 75798, 123648,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 185103, 0, 3,
                                                                       175086, 117168, 175842,
                                                                       75798, 76248, 124323,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 186048, 0, 3,
                                                                       176598, 118248, 177543,
                                                                       77148, 77698, 124998,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 187203, 0, 3,
                                                                       177543, 118923, 178488,
                                                                       77698, 78248, 125823,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 188358, 0, 3,
                                                                       178488, 119598, 179433,
                                                                       78248, 78798, 126648,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 189513, 0, 3,
                                                                       179433, 120273, 180378,
                                                                       78798, 79348, 127473,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 190668, 0, 3,
                                                                       181323, 121623, 182268,
                                                                       80448, 80998, 128298,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 191823, 0, 3,
                                                                       182268, 122298, 183213,
                                                                       80998, 81548, 129123,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 192978, 0, 3,
                                                                       183213, 122973, 184158,
                                                                       81548, 82098, 129948,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 194133, 0, 3,
                                                                       184158, 123648, 185103,
                                                                       82098, 82648, 130773,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 195288, 0, 3,
                                                                       186048, 124998, 187203,
                                                                       83748, 84408, 131598,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 196674, 0, 3,
                                                                       187203, 125823, 188358,
                                                                       84408, 85068, 132588,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 198060, 0, 3,
                                                                       188358, 126648, 189513,
                                                                       85068, 85728, 133578,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 199446, 0, 3,
                                                                       190668, 128298, 191823,
                                                                       87048, 87708, 134568,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 200832, 0, 3,
                                                                       191823, 129123, 192978,
                                                                       87708, 88368, 135558,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 202218, 0, 3,
                                                                       192978, 129948, 194133,
                                                                       88368, 89028, 136548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203604, 3, 90348,
                                                                       90363, 137580, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203632, 3, 90363,
                                                                       90378, 137601, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203660, 3, 90378,
                                                                       90393, 137622, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203688, 3, 90393,
                                                                       90408, 137643, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203716, 3, 90408,
                                                                       90423, 137664, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203744, 3, 90423,
                                                                       90438, 137685, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203772, 3, 90438,
                                                                       90453, 137706, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203800, 3, 90453,
                                                                       90468, 137727, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203828, 3, 90468,
                                                                       90483, 137748, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203856, 3, 90483,
                                                                       90498, 137769, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203884, 3, 90498,
                                                                       90513, 137790, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203912, 3, 90543,
                                                                       90558, 137853, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203940, 3, 90558,
                                                                       90573, 137874, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203968, 3, 90573,
                                                                       90588, 137895, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 203996, 3, 90588,
                                                                       90603, 137916, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204024, 3, 90603,
                                                                       90618, 137937, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204052, 3, 90618,
                                                                       90633, 137958, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204080, 3, 90633,
                                                                       90648, 137979, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204108, 3, 90648,
                                                                       90663, 138000, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204136, 3, 90663,
                                                                       90678, 138021, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204164, 3, 90678,
                                                                       90693, 138042, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 204192, 3, 90693,
                                                                       90708, 138063, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 204220, 0, 3,
                                                                       203604, 137580, 203632,
                                                                       90738, 90783, 138210,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 204304, 0, 3,
                                                                       203632, 137601, 203660,
                                                                       90783, 90828, 138273,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 204388, 0, 3,
                                                                       203660, 137622, 203688,
                                                                       90828, 90873, 138336,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 204472, 0, 3,
                                                                       203688, 137643, 203716,
                                                                       90873, 90918, 138399,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 204556, 0, 3,
                                                                       203716, 137664, 203744,
                                                                       90918, 90963, 138462,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 204640, 0, 3,
                                                                       203744, 137685, 203772,
                                                                       90963, 91008, 138525,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 204724, 0, 3,
                                                                       203772, 137706, 203800,
                                                                       91008, 91053, 138588,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 204808, 0, 3,
                                                                       203800, 137727, 203828,
                                                                       91053, 91098, 138651,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 204892, 0, 3,
                                                                       203828, 137748, 203856,
                                                                       91098, 91143, 138714,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 204976, 0, 3,
                                                                       203856, 137769, 203884,
                                                                       91143, 91188, 138777,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 205060, 0, 3,
                                                                       203912, 137853, 203940,
                                                                       91278, 91323, 138966,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 205144, 0, 3,
                                                                       203940, 137874, 203968,
                                                                       91323, 91368, 139029,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 205228, 0, 3,
                                                                       203968, 137895, 203996,
                                                                       91368, 91413, 139092,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 205312, 0, 3,
                                                                       203996, 137916, 204024,
                                                                       91413, 91458, 139155,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 205396, 0, 3,
                                                                       204024, 137937, 204052,
                                                                       91458, 91503, 139218,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 205480, 0, 3,
                                                                       204052, 137958, 204080,
                                                                       91503, 91548, 139281,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 205564, 0, 3,
                                                                       204080, 137979, 204108,
                                                                       91548, 91593, 139344,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 205648, 0, 3,
                                                                       204108, 138000, 204136,
                                                                       91593, 91638, 139407,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 205732, 0, 3,
                                                                       204136, 138021, 204164,
                                                                       91638, 91683, 139470,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 205816, 0, 3,
                                                                       204164, 138042, 204192,
                                                                       91683, 91728, 139533,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 205900, 0, 3,
                                                                       204220, 138210, 204304,
                                                                       91818, 91908, 139848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 206068, 0, 3,
                                                                       204304, 138273, 204388,
                                                                       91908, 91998, 139974,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 206236, 0, 3,
                                                                       204388, 138336, 204472,
                                                                       91998, 92088, 140100,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 206404, 0, 3,
                                                                       204472, 138399, 204556,
                                                                       92088, 92178, 140226,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 206572, 0, 3,
                                                                       204556, 138462, 204640,
                                                                       92178, 92268, 140352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 206740, 0, 3,
                                                                       204640, 138525, 204724,
                                                                       92268, 92358, 140478,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 206908, 0, 3,
                                                                       204724, 138588, 204808,
                                                                       92358, 92448, 140604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 207076, 0, 3,
                                                                       204808, 138651, 204892,
                                                                       92448, 92538, 140730,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 207244, 0, 3,
                                                                       204892, 138714, 204976,
                                                                       92538, 92628, 140856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 207412, 0, 3,
                                                                       205060, 138966, 205144,
                                                                       92808, 92898, 141234,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 207580, 0, 3,
                                                                       205144, 139029, 205228,
                                                                       92898, 92988, 141360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 207748, 0, 3,
                                                                       205228, 139092, 205312,
                                                                       92988, 93078, 141486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 207916, 0, 3,
                                                                       205312, 139155, 205396,
                                                                       93078, 93168, 141612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 208084, 0, 3,
                                                                       205396, 139218, 205480,
                                                                       93168, 93258, 141738,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 208252, 0, 3,
                                                                       205480, 139281, 205564,
                                                                       93258, 93348, 141864,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 208420, 0, 3,
                                                                       205564, 139344, 205648,
                                                                       93348, 93438, 141990,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 208588, 0, 3,
                                                                       205648, 139407, 205732,
                                                                       93438, 93528, 142116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 208756, 0, 3,
                                                                       205732, 139470, 205816,
                                                                       93528, 93618, 142242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 208924, 0, 3,
                                                                       205900, 139848, 206068,
                                                                       93798, 93948, 142788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 209204, 0, 3,
                                                                       206068, 139974, 206236,
                                                                       93948, 94098, 142998,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 209484, 0, 3,
                                                                       206236, 140100, 206404,
                                                                       94098, 94248, 143208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 209764, 0, 3,
                                                                       206404, 140226, 206572,
                                                                       94248, 94398, 143418,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 210044, 0, 3,
                                                                       206572, 140352, 206740,
                                                                       94398, 94548, 143628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 210324, 0, 3,
                                                                       206740, 140478, 206908,
                                                                       94548, 94698, 143838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 210604, 0, 3,
                                                                       206908, 140604, 207076,
                                                                       94698, 94848, 144048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 210884, 0, 3,
                                                                       207076, 140730, 207244,
                                                                       94848, 94998, 144258,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 211164, 0, 3,
                                                                       207412, 141234, 207580,
                                                                       95298, 95448, 144888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 211444, 0, 3,
                                                                       207580, 141360, 207748,
                                                                       95448, 95598, 145098,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 211724, 0, 3,
                                                                       207748, 141486, 207916,
                                                                       95598, 95748, 145308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 212004, 0, 3,
                                                                       207916, 141612, 208084,
                                                                       95748, 95898, 145518,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 212284, 0, 3,
                                                                       208084, 141738, 208252,
                                                                       95898, 96048, 145728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 212564, 0, 3,
                                                                       208252, 141864, 208420,
                                                                       96048, 96198, 145938,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 212844, 0, 3,
                                                                       208420, 141990, 208588,
                                                                       96198, 96348, 146148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 213124, 0, 3,
                                                                       208588, 142116, 208756,
                                                                       96348, 96498, 146358,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 213404, 0, 3,
                                                                       208924, 142788, 209204,
                                                                       96798, 97023, 147198,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 213824, 0, 3,
                                                                       209204, 142998, 209484,
                                                                       97023, 97248, 147513,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 214244, 0, 3,
                                                                       209484, 143208, 209764,
                                                                       97248, 97473, 147828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 214664, 0, 3,
                                                                       209764, 143418, 210044,
                                                                       97473, 97698, 148143,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 215084, 0, 3,
                                                                       210044, 143628, 210324,
                                                                       97698, 97923, 148458,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 215504, 0, 3,
                                                                       210324, 143838, 210604,
                                                                       97923, 98148, 148773,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 215924, 0, 3,
                                                                       210604, 144048, 210884,
                                                                       98148, 98373, 149088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 216344, 0, 3,
                                                                       211164, 144888, 211444,
                                                                       98823, 99048, 150033,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 216764, 0, 3,
                                                                       211444, 145098, 211724,
                                                                       99048, 99273, 150348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 217184, 0, 3,
                                                                       211724, 145308, 212004,
                                                                       99273, 99498, 150663,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 217604, 0, 3,
                                                                       212004, 145518, 212284,
                                                                       99498, 99723, 150978,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 218024, 0, 3,
                                                                       212284, 145728, 212564,
                                                                       99723, 99948, 151293,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 218444, 0, 3,
                                                                       212564, 145938, 212844,
                                                                       99948, 100173, 151608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 218864, 0, 3,
                                                                       212844, 146148, 213124,
                                                                       100173, 100398, 151923,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 219284, 0, 3,
                                                                       213404, 147198, 213824,
                                                                       100848, 101163, 153120,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 219872, 0, 3,
                                                                       213824, 147513, 214244,
                                                                       101163, 101478, 153561,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 220460, 0, 3,
                                                                       214244, 147828, 214664,
                                                                       101478, 101793, 154002,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 221048, 0, 3,
                                                                       214664, 148143, 215084,
                                                                       101793, 102108, 154443,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 221636, 0, 3,
                                                                       215084, 148458, 215504,
                                                                       102108, 102423, 154884,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 222224, 0, 3,
                                                                       215504, 148773, 215924,
                                                                       102423, 102738, 155325,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 222812, 0, 3,
                                                                       216344, 150033, 216764,
                                                                       103368, 103683, 156648,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 223400, 0, 3,
                                                                       216764, 150348, 217184,
                                                                       103683, 103998, 157089,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 223988, 0, 3,
                                                                       217184, 150663, 217604,
                                                                       103998, 104313, 157530,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 224576, 0, 3,
                                                                       217604, 150978, 218024,
                                                                       104313, 104628, 157971,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 225164, 0, 3,
                                                                       218024, 151293, 218444,
                                                                       104628, 104943, 158412,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 225752, 0, 3,
                                                                       218444, 151608, 218864,
                                                                       104943, 105258, 158853,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 226340, 0, 3,
                                                                       219284, 153120, 219872,
                                                                       105888, 106308, 160470,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 227124, 0, 3,
                                                                       219872, 153561, 220460,
                                                                       106308, 106728, 161058,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 227908, 0, 3,
                                                                       220460, 154002, 221048,
                                                                       106728, 107148, 161646,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 228692, 0, 3,
                                                                       221048, 154443, 221636,
                                                                       107148, 107568, 162234,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 229476, 0, 3,
                                                                       221636, 154884, 222224,
                                                                       107568, 107988, 162822,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 230260, 0, 3,
                                                                       222812, 156648, 223400,
                                                                       108828, 109248, 164586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 231044, 0, 3,
                                                                       223400, 157089, 223988,
                                                                       109248, 109668, 165174,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 231828, 0, 3,
                                                                       223988, 157530, 224576,
                                                                       109668, 110088, 165762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 232612, 0, 3,
                                                                       224576, 157971, 225164,
                                                                       110088, 110508, 166350,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 233396, 0, 3,
                                                                       225164, 158412, 225752,
                                                                       110508, 110928, 166938,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 234180, 0, 3,
                                                                       226340, 160470, 227124,
                                                                       111768, 112308, 169038,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 235188, 0, 3,
                                                                       227124, 161058, 227908,
                                                                       112308, 112848, 169794,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 236196, 0, 3,
                                                                       227908, 161646, 228692,
                                                                       112848, 113388, 170550,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 237204, 0, 3,
                                                                       228692, 162234, 229476,
                                                                       113388, 113928, 171306,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 238212, 0, 3,
                                                                       230260, 164586, 231044,
                                                                       115008, 115548, 173574,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 239220, 0, 3,
                                                                       231044, 165174, 231828,
                                                                       115548, 116088, 174330,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 240228, 0, 3,
                                                                       231828, 165762, 232612,
                                                                       116088, 116628, 175086,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 241236, 0, 3,
                                                                       232612, 166350, 233396,
                                                                       116628, 117168, 175842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 242244, 0, 3,
                                                                       234180, 169038, 235188,
                                                                       118248, 118923, 178488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 243504, 0, 3,
                                                                       235188, 169794, 236196,
                                                                       118923, 119598, 179433,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 244764, 0, 3,
                                                                       236196, 170550, 237204,
                                                                       119598, 120273, 180378,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 246024, 0, 3,
                                                                       238212, 173574, 239220,
                                                                       121623, 122298, 183213,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 247284, 0, 3,
                                                                       239220, 174330, 240228,
                                                                       122298, 122973, 184158,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 248544, 0, 3,
                                                                       240228, 175086, 241236,
                                                                       122973, 123648, 185103,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 249804, 0, 3,
                                                                       242244, 178488, 243504,
                                                                       124998, 125823, 188358,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 251344, 0, 3,
                                                                       243504, 179433, 244764,
                                                                       125823, 126648, 189513,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 252884, 0, 3,
                                                                       246024, 183213, 247284,
                                                                       128298, 129123, 192978,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 254424, 0, 3,
                                                                       247284, 184158, 248544,
                                                                       129123, 129948, 194133,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 255964, 0, 3,
                                                                       249804, 188358, 251344,
                                                                       131598, 132588, 198060,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 257812, 0, 3,
                                                                       252884, 192978, 254424,
                                                                       134568, 135558, 202218,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259660, 3, 137538,
                                                                       137559, 203604, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259696, 3, 137559,
                                                                       137580, 203632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259732, 3, 137580,
                                                                       137601, 203660, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259768, 3, 137601,
                                                                       137622, 203688, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259804, 3, 137622,
                                                                       137643, 203716, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259840, 3, 137643,
                                                                       137664, 203744, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259876, 3, 137664,
                                                                       137685, 203772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259912, 3, 137685,
                                                                       137706, 203800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259948, 3, 137706,
                                                                       137727, 203828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 259984, 3, 137727,
                                                                       137748, 203856, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260020, 3, 137748,
                                                                       137769, 203884, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260056, 3, 137811,
                                                                       137832, 203912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260092, 3, 137832,
                                                                       137853, 203940, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260128, 3, 137853,
                                                                       137874, 203968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260164, 3, 137874,
                                                                       137895, 203996, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260200, 3, 137895,
                                                                       137916, 204024, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260236, 3, 137916,
                                                                       137937, 204052, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260272, 3, 137937,
                                                                       137958, 204080, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260308, 3, 137958,
                                                                       137979, 204108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260344, 3, 137979,
                                                                       138000, 204136, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260380, 3, 138000,
                                                                       138021, 204164, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 260416, 3, 138021,
                                                                       138042, 204192, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 260452, 0, 3,
                                                                       259660, 203604, 259696,
                                                                       138084, 138147, 204220,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 260560, 0, 3,
                                                                       259696, 203632, 259732,
                                                                       138147, 138210, 204304,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 260668, 0, 3,
                                                                       259732, 203660, 259768,
                                                                       138210, 138273, 204388,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 260776, 0, 3,
                                                                       259768, 203688, 259804,
                                                                       138273, 138336, 204472,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 260884, 0, 3,
                                                                       259804, 203716, 259840,
                                                                       138336, 138399, 204556,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 260992, 0, 3,
                                                                       259840, 203744, 259876,
                                                                       138399, 138462, 204640,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 261100, 0, 3,
                                                                       259876, 203772, 259912,
                                                                       138462, 138525, 204724,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 261208, 0, 3,
                                                                       259912, 203800, 259948,
                                                                       138525, 138588, 204808,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 261316, 0, 3,
                                                                       259948, 203828, 259984,
                                                                       138588, 138651, 204892,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 261424, 0, 3,
                                                                       259984, 203856, 260020,
                                                                       138651, 138714, 204976,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 261532, 0, 3,
                                                                       260056, 203912, 260092,
                                                                       138840, 138903, 205060,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 261640, 0, 3,
                                                                       260092, 203940, 260128,
                                                                       138903, 138966, 205144,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 261748, 0, 3,
                                                                       260128, 203968, 260164,
                                                                       138966, 139029, 205228,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 261856, 0, 3,
                                                                       260164, 203996, 260200,
                                                                       139029, 139092, 205312,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 261964, 0, 3,
                                                                       260200, 204024, 260236,
                                                                       139092, 139155, 205396,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 262072, 0, 3,
                                                                       260236, 204052, 260272,
                                                                       139155, 139218, 205480,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 262180, 0, 3,
                                                                       260272, 204080, 260308,
                                                                       139218, 139281, 205564,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 262288, 0, 3,
                                                                       260308, 204108, 260344,
                                                                       139281, 139344, 205648,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 262396, 0, 3,
                                                                       260344, 204136, 260380,
                                                                       139344, 139407, 205732,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 262504, 0, 3,
                                                                       260380, 204164, 260416,
                                                                       139407, 139470, 205816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 262612, 0, 3,
                                                                       260452, 204220, 260560,
                                                                       139596, 139722, 205900,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 262828, 0, 3,
                                                                       260560, 204304, 260668,
                                                                       139722, 139848, 206068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 263044, 0, 3,
                                                                       260668, 204388, 260776,
                                                                       139848, 139974, 206236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 263260, 0, 3,
                                                                       260776, 204472, 260884,
                                                                       139974, 140100, 206404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 263476, 0, 3,
                                                                       260884, 204556, 260992,
                                                                       140100, 140226, 206572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 263692, 0, 3,
                                                                       260992, 204640, 261100,
                                                                       140226, 140352, 206740,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 263908, 0, 3,
                                                                       261100, 204724, 261208,
                                                                       140352, 140478, 206908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 264124, 0, 3,
                                                                       261208, 204808, 261316,
                                                                       140478, 140604, 207076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 264340, 0, 3,
                                                                       261316, 204892, 261424,
                                                                       140604, 140730, 207244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 264556, 0, 3,
                                                                       261532, 205060, 261640,
                                                                       140982, 141108, 207412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 264772, 0, 3,
                                                                       261640, 205144, 261748,
                                                                       141108, 141234, 207580,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 264988, 0, 3,
                                                                       261748, 205228, 261856,
                                                                       141234, 141360, 207748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 265204, 0, 3,
                                                                       261856, 205312, 261964,
                                                                       141360, 141486, 207916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 265420, 0, 3,
                                                                       261964, 205396, 262072,
                                                                       141486, 141612, 208084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 265636, 0, 3,
                                                                       262072, 205480, 262180,
                                                                       141612, 141738, 208252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 265852, 0, 3,
                                                                       262180, 205564, 262288,
                                                                       141738, 141864, 208420,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 266068, 0, 3,
                                                                       262288, 205648, 262396,
                                                                       141864, 141990, 208588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 266284, 0, 3,
                                                                       262396, 205732, 262504,
                                                                       141990, 142116, 208756,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 266500, 0, 3,
                                                                       262612, 205900, 262828,
                                                                       142368, 142578, 208924,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 266860, 0, 3,
                                                                       262828, 206068, 263044,
                                                                       142578, 142788, 209204,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 267220, 0, 3,
                                                                       263044, 206236, 263260,
                                                                       142788, 142998, 209484,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 267580, 0, 3,
                                                                       263260, 206404, 263476,
                                                                       142998, 143208, 209764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 267940, 0, 3,
                                                                       263476, 206572, 263692,
                                                                       143208, 143418, 210044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 268300, 0, 3,
                                                                       263692, 206740, 263908,
                                                                       143418, 143628, 210324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 268660, 0, 3,
                                                                       263908, 206908, 264124,
                                                                       143628, 143838, 210604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 269020, 0, 3,
                                                                       264124, 207076, 264340,
                                                                       143838, 144048, 210884,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 269380, 0, 3,
                                                                       264556, 207412, 264772,
                                                                       144468, 144678, 211164,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 269740, 0, 3,
                                                                       264772, 207580, 264988,
                                                                       144678, 144888, 211444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 270100, 0, 3,
                                                                       264988, 207748, 265204,
                                                                       144888, 145098, 211724,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 270460, 0, 3,
                                                                       265204, 207916, 265420,
                                                                       145098, 145308, 212004,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 270820, 0, 3,
                                                                       265420, 208084, 265636,
                                                                       145308, 145518, 212284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 271180, 0, 3,
                                                                       265636, 208252, 265852,
                                                                       145518, 145728, 212564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 271540, 0, 3,
                                                                       265852, 208420, 266068,
                                                                       145728, 145938, 212844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 271900, 0, 3,
                                                                       266068, 208588, 266284,
                                                                       145938, 146148, 213124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 272260, 0, 3,
                                                                       266500, 208924, 266860,
                                                                       146568, 146883, 213404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 272800, 0, 3,
                                                                       266860, 209204, 267220,
                                                                       146883, 147198, 213824,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 273340, 0, 3,
                                                                       267220, 209484, 267580,
                                                                       147198, 147513, 214244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 273880, 0, 3,
                                                                       267580, 209764, 267940,
                                                                       147513, 147828, 214664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 274420, 0, 3,
                                                                       267940, 210044, 268300,
                                                                       147828, 148143, 215084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 274960, 0, 3,
                                                                       268300, 210324, 268660,
                                                                       148143, 148458, 215504,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 275500, 0, 3,
                                                                       268660, 210604, 269020,
                                                                       148458, 148773, 215924,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 276040, 0, 3,
                                                                       269380, 211164, 269740,
                                                                       149403, 149718, 216344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 276580, 0, 3,
                                                                       269740, 211444, 270100,
                                                                       149718, 150033, 216764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 277120, 0, 3,
                                                                       270100, 211724, 270460,
                                                                       150033, 150348, 217184,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 277660, 0, 3,
                                                                       270460, 212004, 270820,
                                                                       150348, 150663, 217604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 278200, 0, 3,
                                                                       270820, 212284, 271180,
                                                                       150663, 150978, 218024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 278740, 0, 3,
                                                                       271180, 212564, 271540,
                                                                       150978, 151293, 218444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 279280, 0, 3,
                                                                       271540, 212844, 271900,
                                                                       151293, 151608, 218864,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 279820, 0, 3,
                                                                       272260, 213404, 272800,
                                                                       152238, 152679, 219284,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 280576, 0, 3,
                                                                       272800, 213824, 273340,
                                                                       152679, 153120, 219872,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 281332, 0, 3,
                                                                       273340, 214244, 273880,
                                                                       153120, 153561, 220460,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 282088, 0, 3,
                                                                       273880, 214664, 274420,
                                                                       153561, 154002, 221048,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 282844, 0, 3,
                                                                       274420, 215084, 274960,
                                                                       154002, 154443, 221636,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 283600, 0, 3,
                                                                       274960, 215504, 275500,
                                                                       154443, 154884, 222224,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 284356, 0, 3,
                                                                       276040, 216344, 276580,
                                                                       155766, 156207, 222812,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 285112, 0, 3,
                                                                       276580, 216764, 277120,
                                                                       156207, 156648, 223400,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 285868, 0, 3,
                                                                       277120, 217184, 277660,
                                                                       156648, 157089, 223988,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 286624, 0, 3,
                                                                       277660, 217604, 278200,
                                                                       157089, 157530, 224576,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 287380, 0, 3,
                                                                       278200, 218024, 278740,
                                                                       157530, 157971, 225164,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 288136, 0, 3,
                                                                       278740, 218444, 279280,
                                                                       157971, 158412, 225752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 288892, 0, 3,
                                                                       279820, 219284, 280576,
                                                                       159294, 159882, 226340,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 289900, 0, 3,
                                                                       280576, 219872, 281332,
                                                                       159882, 160470, 227124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 290908, 0, 3,
                                                                       281332, 220460, 282088,
                                                                       160470, 161058, 227908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 291916, 0, 3,
                                                                       282088, 221048, 282844,
                                                                       161058, 161646, 228692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 292924, 0, 3,
                                                                       282844, 221636, 283600,
                                                                       161646, 162234, 229476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 293932, 0, 3,
                                                                       284356, 222812, 285112,
                                                                       163410, 163998, 230260,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 294940, 0, 3,
                                                                       285112, 223400, 285868,
                                                                       163998, 164586, 231044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 295948, 0, 3,
                                                                       285868, 223988, 286624,
                                                                       164586, 165174, 231828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 296956, 0, 3,
                                                                       286624, 224576, 287380,
                                                                       165174, 165762, 232612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 297964, 0, 3,
                                                                       287380, 225164, 288136,
                                                                       165762, 166350, 233396,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 298972, 0, 3,
                                                                       288892, 226340, 289900,
                                                                       167526, 168282, 234180,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 300268, 0, 3,
                                                                       289900, 227124, 290908,
                                                                       168282, 169038, 235188,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 301564, 0, 3,
                                                                       290908, 227908, 291916,
                                                                       169038, 169794, 236196,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 302860, 0, 3,
                                                                       291916, 228692, 292924,
                                                                       169794, 170550, 237204,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 304156, 0, 3,
                                                                       293932, 230260, 294940,
                                                                       172062, 172818, 238212,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 305452, 0, 3,
                                                                       294940, 231044, 295948,
                                                                       172818, 173574, 239220,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 306748, 0, 3,
                                                                       295948, 231828, 296956,
                                                                       173574, 174330, 240228,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 308044, 0, 3,
                                                                       296956, 232612, 297964,
                                                                       174330, 175086, 241236,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 309340, 0, 3,
                                                                       298972, 234180, 300268,
                                                                       176598, 177543, 242244,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 310960, 0, 3,
                                                                       300268, 235188, 301564,
                                                                       177543, 178488, 243504,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 312580, 0, 3,
                                                                       301564, 236196, 302860,
                                                                       178488, 179433, 244764,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 314200, 0, 3,
                                                                       304156, 238212, 305452,
                                                                       181323, 182268, 246024,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 315820, 0, 3,
                                                                       305452, 239220, 306748,
                                                                       182268, 183213, 247284,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 317440, 0, 3,
                                                                       306748, 240228, 308044,
                                                                       183213, 184158, 248544,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 319060, 0, 3,
                                                                       309340, 242244, 310960,
                                                                       186048, 187203, 249804,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 321040, 0, 3,
                                                                       310960, 243504, 312580,
                                                                       187203, 188358, 251344,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 323020, 0, 3,
                                                                       314200, 246024, 315820,
                                                                       190668, 191823, 252884,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 325000, 0, 3,
                                                                       315820, 247284, 317440,
                                                                       191823, 192978, 254424,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 326980, 0, 3,
                                                                       319060, 249804, 321040,
                                                                       195288, 196674, 255964,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 329356, 0, 3,
                                                                       323020, 252884, 325000,
                                                                       199446, 200832, 257812,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 331732, 288892, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 333160, 293932, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 334588, 298972, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 336424, 304156, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 338260, 309340, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 340555, 314200, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 342850, 319060, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 345655, 323020, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 348460, 326980, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 351826, 329356, 2376, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 332740, 331732, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 334168, 333160, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 335884, 334588, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 337720, 336424, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 339880, 338260, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 342175, 340555, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 344830, 342850, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 347635, 345655, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 350836, 348460, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 354202, 351826, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 355192, 332740, 335884, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 356452, 334168, 337720, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 357712, 335884, 339880, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 359332, 337720, 342175, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 360952, 339880, 344830, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 362977, 342175, 347635, 15, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 365002, 344830, 350836, 15, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 367477, 347635, 354202, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 369952, 355192, 357712, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 372472, 356452, 359332, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 374992, 357712, 360952, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 378232, 359332, 362977, 15, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 381472, 360952, 365002, 15, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 385522, 362977, 367477, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 389572, 369952, 374992, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 393772, 372472, 378232, 15, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 397972, 374992, 381472, 15, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 403372, 378232, 385522, 15, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 408772, 389572, 397972, 15, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 415072, 393772, 403372, 15, nmax);

        simdtrf::transform_i_inner(buffer, 421372, 415072, 15, 15, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 421372, 195, nmax);

        simdtrf::transform_i_inner(buffer, 421372, 408772, 15, 15, nmax);

        simdtrf::transform_g_outer(values + 1755 * nvalues + n * npairs, nvalues, buffer, 421372,
                                   195, nmax);
    }

    for (size_t m = 0; m < 3510; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
