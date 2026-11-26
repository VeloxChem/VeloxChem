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


#include "SimdThreeCenterElectronRepulsionRsRecHIF.hpp"

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
compute_rs_hif_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hif_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 133147, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2002 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 133147, 60396, 9926, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                                                            14}, ncols, fj, i * nprim_b + j, fq,
                                                            omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 21, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
                                                        ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 69, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 72, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 75, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 78, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 81, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 84, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 87, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 90, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 93, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 96, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 99, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 102, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 105, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 108, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 111, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 114, 0, 3, 7, 8,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 120, 0, 3, 8, 9,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 126, 0, 3, 9, 10,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 132, 0, 3, 10, 11,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 138, 0, 3, 11, 12,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 144, 0, 3, 12, 13,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 150, 0, 3, 13, 14,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 156, 0, 3, 14, 15,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 162, 0, 3, 15, 16,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 168, 0, 3, 16, 17,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 174, 0, 3, 17, 18,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 180, 0, 3, 18, 19,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 186, 0, 3, 22, 23,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 192, 0, 3, 23, 24,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 198, 0, 3, 24, 25,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 204, 0, 3, 25, 26,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 210, 0, 3, 26, 27,
                                                                       87, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 216, 0, 3, 27, 28,
                                                                       90, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 222, 0, 3, 28, 29,
                                                                       93, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 228, 0, 3, 29, 30,
                                                                       96, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 234, 0, 3, 30, 31,
                                                                       99, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 240, 0, 3, 31, 32,
                                                                       102, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 246, 0, 3, 32, 33,
                                                                       105, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 252, 0, 3, 33, 34,
                                                                       108, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 36, 39,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 39, 42,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 42, 45,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 45, 48,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 48, 51,
                                                                       138, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 51, 54,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 54, 57,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 57, 60,
                                                                       156, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 60, 63,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 63, 66,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 66, 69,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 75, 78,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 78, 81,
                                                                       192, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 81, 84,
                                                                       198, 204, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 84, 87,
                                                                       204, 210, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 87, 90,
                                                                       210, 216, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 90, 93,
                                                                       216, 222, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 93, 96,
                                                                       222, 228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 96, 99,
                                                                       228, 234, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 99,
                                                                       102, 234, 240, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 102,
                                                                       105, 240, 246, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 105,
                                                                       108, 246, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 114,
                                                                       120, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 493, 0, 3, 120,
                                                                       126, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 126,
                                                                       132, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 523, 0, 3, 132,
                                                                       138, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 538, 0, 3, 138,
                                                                       144, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 553, 0, 3, 144,
                                                                       150, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 568, 0, 3, 150,
                                                                       156, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 583, 0, 3, 156,
                                                                       162, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 598, 0, 3, 162,
                                                                       168, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 613, 0, 3, 168,
                                                                       174, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 628, 0, 3, 186,
                                                                       192, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 643, 0, 3, 192,
                                                                       198, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 658, 0, 3, 198,
                                                                       204, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 673, 0, 3, 204,
                                                                       210, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 688, 0, 3, 210,
                                                                       216, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 703, 0, 3, 216,
                                                                       222, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 718, 0, 3, 222,
                                                                       228, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 733, 0, 3, 228,
                                                                       234, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 748, 0, 3, 234,
                                                                       240, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 763, 0, 3, 240,
                                                                       246, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 778, 0, 3, 258,
                                                                       268, 478, 493, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 799, 0, 3, 268,
                                                                       278, 493, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 820, 0, 3, 278,
                                                                       288, 508, 523, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 841, 0, 3, 288,
                                                                       298, 523, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 862, 0, 3, 298,
                                                                       308, 538, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 883, 0, 3, 308,
                                                                       318, 553, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 904, 0, 3, 318,
                                                                       328, 568, 583, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 925, 0, 3, 328,
                                                                       338, 583, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 946, 0, 3, 338,
                                                                       348, 598, 613, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 967, 0, 3, 368,
                                                                       378, 628, 643, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 988, 0, 3, 378,
                                                                       388, 643, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1009, 0, 3, 388,
                                                                       398, 658, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 398,
                                                                       408, 673, 688, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1051, 0, 3, 408,
                                                                       418, 688, 703, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 418,
                                                                       428, 703, 718, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1093, 0, 3, 428,
                                                                       438, 718, 733, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 438,
                                                                       448, 733, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1135, 0, 3, 448,
                                                                       458, 748, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 478,
                                                                       493, 778, 799, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 493,
                                                                       508, 799, 820, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 508,
                                                                       523, 820, 841, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 523,
                                                                       538, 841, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 538,
                                                                       553, 862, 883, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 553,
                                                                       568, 883, 904, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 568,
                                                                       583, 904, 925, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 583,
                                                                       598, 925, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 628,
                                                                       643, 967, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 643,
                                                                       658, 988, 1009, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 658,
                                                                       673, 1009, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 673,
                                                                       688, 1030, 1051, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 688,
                                                                       703, 1051, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 703,
                                                                       718, 1072, 1093, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 718,
                                                                       733, 1093, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 733,
                                                                       748, 1114, 1135, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 778,
                                                                       799, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1640, 0, 3, 799,
                                                                       820, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1676, 0, 3, 820,
                                                                       841, 1212, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1712, 0, 3, 841,
                                                                       862, 1240, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1748, 0, 3, 862,
                                                                       883, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1784, 0, 3, 883,
                                                                       904, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1820, 0, 3, 904,
                                                                       925, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 967,
                                                                       988, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1892, 0, 3, 988,
                                                                       1009, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1928, 0, 3, 1009,
                                                                       1030, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1964, 0, 3, 1030,
                                                                       1051, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2000, 0, 3, 1051,
                                                                       1072, 1492, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2036, 0, 3, 1072,
                                                                       1093, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2072, 0, 3, 1093,
                                                                       1114, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 1156,
                                                                       1184, 1604, 1640, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1184,
                                                                       1212, 1640, 1676, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2198, 0, 3, 1212,
                                                                       1240, 1676, 1712, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2243, 0, 3, 1240,
                                                                       1268, 1712, 1748, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2288, 0, 3, 1268,
                                                                       1296, 1748, 1784, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2333, 0, 3, 1296,
                                                                       1324, 1784, 1820, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2378, 0, 3, 1380,
                                                                       1408, 1856, 1892, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2423, 0, 3, 1408,
                                                                       1436, 1892, 1928, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2468, 0, 3, 1436,
                                                                       1464, 1928, 1964, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2513, 0, 3, 1464,
                                                                       1492, 1964, 2000, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2558, 0, 3, 1492,
                                                                       1520, 2000, 2036, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2603, 0, 3, 1520,
                                                                       1548, 2036, 2072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1604,
                                                                       1640, 2108, 2153, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2703, 0, 3, 1640,
                                                                       1676, 2153, 2198, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2758, 0, 3, 1676,
                                                                       1712, 2198, 2243, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2813, 0, 3, 1712,
                                                                       1748, 2243, 2288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2868, 0, 3, 1748,
                                                                       1784, 2288, 2333, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2923, 0, 3, 1856,
                                                                       1892, 2378, 2423, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2978, 0, 3, 1892,
                                                                       1928, 2423, 2468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3033, 0, 3, 1928,
                                                                       1964, 2468, 2513, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 1964,
                                                                       2000, 2513, 2558, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3143, 0, 3, 2000,
                                                                       2036, 2558, 2603, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3198, 0, 3, 2108,
                                                                       2153, 2648, 2703, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3264, 0, 3, 2153,
                                                                       2198, 2703, 2758, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3330, 0, 3, 2198,
                                                                       2243, 2758, 2813, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3396, 0, 3, 2243,
                                                                       2288, 2813, 2868, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3462, 0, 3, 2378,
                                                                       2423, 2923, 2978, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3528, 0, 3, 2423,
                                                                       2468, 2978, 3033, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3594, 0, 3, 2468,
                                                                       2513, 3033, 3088, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3660, 0, 3, 2513,
                                                                       2558, 3088, 3143, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3726, 0, 3, 2648,
                                                                       2703, 3198, 3264, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3804, 0, 3, 2703,
                                                                       2758, 3264, 3330, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3882, 0, 3, 2758,
                                                                       2813, 3330, 3396, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3960, 0, 3, 2923,
                                                                       2978, 3462, 3528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 4038, 0, 3, 2978,
                                                                       3033, 3528, 3594, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 4116, 0, 3, 3033,
                                                                       3088, 3594, 3660, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4194, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4197, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4200, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4203, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4206, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4209, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4212, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4215, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4218, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4221, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4224, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4227, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4230, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4233, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4236, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4239, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4242, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4245, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4248, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4251, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4254, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4257, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4260, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4263, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4266, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4269, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4272, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4275, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4278, 3, 9, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4287, 3, 10, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4296, 3, 11, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4305, 3, 12, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4314, 3, 13, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4323, 3, 14, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4332, 3, 15, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4341, 3, 16, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4350, 3, 17, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4359, 3, 18, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4368, 3, 19, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4377, 3, 24, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4386, 3, 25, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4395, 3, 26, 87,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4404, 3, 27, 90,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4413, 3, 28, 93,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4422, 3, 29, 96,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4431, 3, 30, 99,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4440, 3, 31, 102,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4449, 3, 32, 105,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4458, 3, 33, 108,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4467, 3, 34, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4476, 3, 36, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4494, 3, 39, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4512, 3, 42, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4530, 3, 45, 132,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4548, 3, 48, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4566, 3, 51, 144,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4584, 3, 54, 150,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4602, 3, 57, 156,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4620, 3, 60, 162,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4638, 3, 63, 168,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4656, 3, 66, 174,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4674, 3, 69, 180,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4692, 3, 75, 186,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4710, 3, 78, 192,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4728, 3, 81, 198,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4746, 3, 84, 204,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4764, 3, 87, 210,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4782, 3, 90, 216,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4800, 3, 93, 222,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4818, 3, 96, 228,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4836, 3, 99, 234,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4854, 3, 102, 240,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4872, 3, 105, 246,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4890, 3, 108, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4908, 3, 114, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4938, 3, 120, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4968, 3, 126, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4998, 3, 132, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5028, 3, 138, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5058, 3, 144, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5088, 3, 150, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5118, 3, 156, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5148, 3, 162, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5178, 3, 168, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5208, 3, 174, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5238, 3, 186, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5268, 3, 192, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5298, 3, 198, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5328, 3, 204, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5358, 3, 210, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5388, 3, 216, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5418, 3, 222, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5448, 3, 228, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5478, 3, 234, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5508, 3, 240, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5538, 3, 246, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5568, 3, 258, 478,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5613, 3, 268, 493,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5658, 3, 278, 508,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5703, 3, 288, 523,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5748, 3, 298, 538,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5793, 3, 308, 553,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5838, 3, 318, 568,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5883, 3, 328, 583,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5928, 3, 338, 598,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5973, 3, 348, 613,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6018, 3, 368, 628,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6063, 3, 378, 643,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6108, 3, 388, 658,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6153, 3, 398, 673,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6198, 3, 408, 688,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6243, 3, 418, 703,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6288, 3, 428, 718,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6333, 3, 438, 733,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6378, 3, 448, 748,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6423, 3, 458, 763,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6468, 3, 478, 778,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6531, 3, 493, 799,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6594, 3, 508, 820,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6657, 3, 523, 841,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6720, 3, 538, 862,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6783, 3, 553, 883,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6846, 3, 568, 904,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6909, 3, 583, 925,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6972, 3, 598, 946,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7035, 3, 628, 967,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7098, 3, 643, 988,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7161, 3, 658,
                                                                       1009, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7224, 3, 673,
                                                                       1030, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7287, 3, 688,
                                                                       1051, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7350, 3, 703,
                                                                       1072, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7413, 3, 718,
                                                                       1093, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7476, 3, 733,
                                                                       1114, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7539, 3, 748,
                                                                       1135, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7602, 3, 778,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7686, 3, 799,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7770, 3, 820,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7854, 3, 841,
                                                                       1240, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7938, 3, 862,
                                                                       1268, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8022, 3, 883,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8106, 3, 904,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8190, 3, 925,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8274, 3, 967,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8358, 3, 988,
                                                                       1408, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8442, 3, 1009,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8526, 3, 1030,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8610, 3, 1051,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8694, 3, 1072,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8778, 3, 1093,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8862, 3, 1114,
                                                                       1576, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8946, 3, 1156,
                                                                       1604, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9054, 3, 1184,
                                                                       1640, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9162, 3, 1212,
                                                                       1676, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9270, 3, 1240,
                                                                       1712, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9378, 3, 1268,
                                                                       1748, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9486, 3, 1296,
                                                                       1784, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9594, 3, 1324,
                                                                       1820, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9702, 3, 1380,
                                                                       1856, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9810, 3, 1408,
                                                                       1892, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9918, 3, 1436,
                                                                       1928, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10026, 3, 1464,
                                                                       1964, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10134, 3, 1492,
                                                                       2000, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10242, 3, 1520,
                                                                       2036, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10350, 3, 1548,
                                                                       2072, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10458, 3, 1604,
                                                                       2108, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10593, 3, 1640,
                                                                       2153, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10728, 3, 1676,
                                                                       2198, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10863, 3, 1712,
                                                                       2243, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10998, 3, 1748,
                                                                       2288, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11133, 3, 1784,
                                                                       2333, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11268, 3, 1856,
                                                                       2378, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11403, 3, 1892,
                                                                       2423, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11538, 3, 1928,
                                                                       2468, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11673, 3, 1964,
                                                                       2513, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11808, 3, 2000,
                                                                       2558, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11943, 3, 2036,
                                                                       2603, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12078, 3, 2108,
                                                                       2648, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12243, 3, 2153,
                                                                       2703, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12408, 3, 2198,
                                                                       2758, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12573, 3, 2243,
                                                                       2813, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12738, 3, 2288,
                                                                       2868, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12903, 3, 2378,
                                                                       2923, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13068, 3, 2423,
                                                                       2978, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13233, 3, 2468,
                                                                       3033, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13398, 3, 2513,
                                                                       3088, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13563, 3, 2558,
                                                                       3143, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 13728, 3, 2648,
                                                                       3198, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 13926, 3, 2703,
                                                                       3264, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 14124, 3, 2758,
                                                                       3330, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 14322, 3, 2813,
                                                                       3396, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 14520, 3, 2923,
                                                                       3462, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 14718, 3, 2978,
                                                                       3528, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 14916, 3, 3033,
                                                                       3594, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 15114, 3, 3088,
                                                                       3660, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 15312, 3, 3198,
                                                                       3726, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 15546, 3, 3264,
                                                                       3804, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 15780, 3, 3330,
                                                                       3882, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 16014, 3, 3462,
                                                                       3960, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 16248, 3, 3528,
                                                                       4038, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 16482, 3, 3594,
                                                                       4116, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16716, 3, 7, 8,
                                                                       4200, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16722, 3, 8, 9,
                                                                       4203, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16728, 3, 9, 10,
                                                                       4206, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16734, 3, 10, 11,
                                                                       4209, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16740, 3, 11, 12,
                                                                       4212, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16746, 3, 12, 13,
                                                                       4215, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16752, 3, 13, 14,
                                                                       4218, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16758, 3, 14, 15,
                                                                       4221, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16764, 3, 15, 16,
                                                                       4224, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16770, 3, 16, 17,
                                                                       4227, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16776, 3, 17, 18,
                                                                       4230, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16782, 3, 18, 19,
                                                                       4233, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16788, 3, 22, 23,
                                                                       4242, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16794, 3, 23, 24,
                                                                       4245, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16800, 3, 24, 25,
                                                                       4248, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16806, 3, 25, 26,
                                                                       4251, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16812, 3, 26, 27,
                                                                       4254, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16818, 3, 27, 28,
                                                                       4257, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16824, 3, 28, 29,
                                                                       4260, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16830, 3, 29, 30,
                                                                       4263, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16836, 3, 30, 31,
                                                                       4266, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16842, 3, 31, 32,
                                                                       4269, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16848, 3, 32, 33,
                                                                       4272, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16854, 3, 33, 34,
                                                                       4275, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16860, 0, 3,
                                                                       16716, 4200, 16722, 4278,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16878, 0, 3,
                                                                       16722, 4203, 16728, 4287,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16896, 0, 3,
                                                                       16728, 4206, 16734, 4296,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16914, 0, 3,
                                                                       16734, 4209, 16740, 4305,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16932, 0, 3,
                                                                       16740, 4212, 16746, 4314,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16950, 0, 3,
                                                                       16746, 4215, 16752, 4323,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16968, 0, 3,
                                                                       16752, 4218, 16758, 4332,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16986, 0, 3,
                                                                       16758, 4221, 16764, 4341,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17004, 0, 3,
                                                                       16764, 4224, 16770, 4350,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17022, 0, 3,
                                                                       16770, 4227, 16776, 4359,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17040, 0, 3,
                                                                       16776, 4230, 16782, 4368,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17058, 0, 3,
                                                                       16788, 4242, 16794, 4377,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17076, 0, 3,
                                                                       16794, 4245, 16800, 4386,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17094, 0, 3,
                                                                       16800, 4248, 16806, 4395,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17112, 0, 3,
                                                                       16806, 4251, 16812, 4404,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17130, 0, 3,
                                                                       16812, 4254, 16818, 4413,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17148, 0, 3,
                                                                       16818, 4257, 16824, 4422,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17166, 0, 3,
                                                                       16824, 4260, 16830, 4431,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17184, 0, 3,
                                                                       16830, 4263, 16836, 4440,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17202, 0, 3,
                                                                       16836, 4266, 16842, 4449,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17220, 0, 3,
                                                                       16842, 4269, 16848, 4458,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 17238, 0, 3,
                                                                       16848, 4272, 16854, 4467,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17256, 0, 3,
                                                                       16860, 4278, 16878, 114,
                                                                       120, 4512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17292, 0, 3,
                                                                       16878, 4287, 16896, 120,
                                                                       126, 4530, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17328, 0, 3,
                                                                       16896, 4296, 16914, 126,
                                                                       132, 4548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17364, 0, 3,
                                                                       16914, 4305, 16932, 132,
                                                                       138, 4566, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17400, 0, 3,
                                                                       16932, 4314, 16950, 138,
                                                                       144, 4584, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17436, 0, 3,
                                                                       16950, 4323, 16968, 144,
                                                                       150, 4602, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17472, 0, 3,
                                                                       16968, 4332, 16986, 150,
                                                                       156, 4620, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17508, 0, 3,
                                                                       16986, 4341, 17004, 156,
                                                                       162, 4638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17544, 0, 3,
                                                                       17004, 4350, 17022, 162,
                                                                       168, 4656, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17580, 0, 3,
                                                                       17022, 4359, 17040, 168,
                                                                       174, 4674, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17616, 0, 3,
                                                                       17058, 4377, 17076, 186,
                                                                       192, 4728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17652, 0, 3,
                                                                       17076, 4386, 17094, 192,
                                                                       198, 4746, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17688, 0, 3,
                                                                       17094, 4395, 17112, 198,
                                                                       204, 4764, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17724, 0, 3,
                                                                       17112, 4404, 17130, 204,
                                                                       210, 4782, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17760, 0, 3,
                                                                       17130, 4413, 17148, 210,
                                                                       216, 4800, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17796, 0, 3,
                                                                       17148, 4422, 17166, 216,
                                                                       222, 4818, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17832, 0, 3,
                                                                       17166, 4431, 17184, 222,
                                                                       228, 4836, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17868, 0, 3,
                                                                       17184, 4440, 17202, 228,
                                                                       234, 4854, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17904, 0, 3,
                                                                       17202, 4449, 17220, 234,
                                                                       240, 4872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17940, 0, 3,
                                                                       17220, 4458, 17238, 240,
                                                                       246, 4890, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 17976, 0, 3,
                                                                       17256, 4512, 17292, 258,
                                                                       268, 4968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18036, 0, 3,
                                                                       17292, 4530, 17328, 268,
                                                                       278, 4998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18096, 0, 3,
                                                                       17328, 4548, 17364, 278,
                                                                       288, 5028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18156, 0, 3,
                                                                       17364, 4566, 17400, 288,
                                                                       298, 5058, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18216, 0, 3,
                                                                       17400, 4584, 17436, 298,
                                                                       308, 5088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18276, 0, 3,
                                                                       17436, 4602, 17472, 308,
                                                                       318, 5118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18336, 0, 3,
                                                                       17472, 4620, 17508, 318,
                                                                       328, 5148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18396, 0, 3,
                                                                       17508, 4638, 17544, 328,
                                                                       338, 5178, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18456, 0, 3,
                                                                       17544, 4656, 17580, 338,
                                                                       348, 5208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18516, 0, 3,
                                                                       17616, 4728, 17652, 368,
                                                                       378, 5298, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18576, 0, 3,
                                                                       17652, 4746, 17688, 378,
                                                                       388, 5328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18636, 0, 3,
                                                                       17688, 4764, 17724, 388,
                                                                       398, 5358, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18696, 0, 3,
                                                                       17724, 4782, 17760, 398,
                                                                       408, 5388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18756, 0, 3,
                                                                       17760, 4800, 17796, 408,
                                                                       418, 5418, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18816, 0, 3,
                                                                       17796, 4818, 17832, 418,
                                                                       428, 5448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18876, 0, 3,
                                                                       17832, 4836, 17868, 428,
                                                                       438, 5478, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18936, 0, 3,
                                                                       17868, 4854, 17904, 438,
                                                                       448, 5508, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18996, 0, 3,
                                                                       17904, 4872, 17940, 448,
                                                                       458, 5538, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19056, 0, 3,
                                                                       17976, 4968, 18036, 478,
                                                                       493, 5658, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19146, 0, 3,
                                                                       18036, 4998, 18096, 493,
                                                                       508, 5703, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19236, 0, 3,
                                                                       18096, 5028, 18156, 508,
                                                                       523, 5748, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19326, 0, 3,
                                                                       18156, 5058, 18216, 523,
                                                                       538, 5793, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19416, 0, 3,
                                                                       18216, 5088, 18276, 538,
                                                                       553, 5838, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19506, 0, 3,
                                                                       18276, 5118, 18336, 553,
                                                                       568, 5883, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19596, 0, 3,
                                                                       18336, 5148, 18396, 568,
                                                                       583, 5928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19686, 0, 3,
                                                                       18396, 5178, 18456, 583,
                                                                       598, 5973, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19776, 0, 3,
                                                                       18516, 5298, 18576, 628,
                                                                       643, 6108, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19866, 0, 3,
                                                                       18576, 5328, 18636, 643,
                                                                       658, 6153, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19956, 0, 3,
                                                                       18636, 5358, 18696, 658,
                                                                       673, 6198, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20046, 0, 3,
                                                                       18696, 5388, 18756, 673,
                                                                       688, 6243, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20136, 0, 3,
                                                                       18756, 5418, 18816, 688,
                                                                       703, 6288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20226, 0, 3,
                                                                       18816, 5448, 18876, 703,
                                                                       718, 6333, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20316, 0, 3,
                                                                       18876, 5478, 18936, 718,
                                                                       733, 6378, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20406, 0, 3,
                                                                       18936, 5508, 18996, 733,
                                                                       748, 6423, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 20496, 0, 3,
                                                                       19056, 5658, 19146, 778,
                                                                       799, 6594, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 20622, 0, 3,
                                                                       19146, 5703, 19236, 799,
                                                                       820, 6657, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 20748, 0, 3,
                                                                       19236, 5748, 19326, 820,
                                                                       841, 6720, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 20874, 0, 3,
                                                                       19326, 5793, 19416, 841,
                                                                       862, 6783, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21000, 0, 3,
                                                                       19416, 5838, 19506, 862,
                                                                       883, 6846, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21126, 0, 3,
                                                                       19506, 5883, 19596, 883,
                                                                       904, 6909, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21252, 0, 3,
                                                                       19596, 5928, 19686, 904,
                                                                       925, 6972, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21378, 0, 3,
                                                                       19776, 6108, 19866, 967,
                                                                       988, 7161, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21504, 0, 3,
                                                                       19866, 6153, 19956, 988,
                                                                       1009, 7224, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21630, 0, 3,
                                                                       19956, 6198, 20046, 1009,
                                                                       1030, 7287, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21756, 0, 3,
                                                                       20046, 6243, 20136, 1030,
                                                                       1051, 7350, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21882, 0, 3,
                                                                       20136, 6288, 20226, 1051,
                                                                       1072, 7413, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22008, 0, 3,
                                                                       20226, 6333, 20316, 1072,
                                                                       1093, 7476, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22134, 0, 3,
                                                                       20316, 6378, 20406, 1093,
                                                                       1114, 7539, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22260, 0, 3,
                                                                       20496, 6594, 20622, 1156,
                                                                       1184, 7770, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22428, 0, 3,
                                                                       20622, 6657, 20748, 1184,
                                                                       1212, 7854, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22596, 0, 3,
                                                                       20748, 6720, 20874, 1212,
                                                                       1240, 7938, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22764, 0, 3,
                                                                       20874, 6783, 21000, 1240,
                                                                       1268, 8022, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 22932, 0, 3,
                                                                       21000, 6846, 21126, 1268,
                                                                       1296, 8106, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23100, 0, 3,
                                                                       21126, 6909, 21252, 1296,
                                                                       1324, 8190, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23268, 0, 3,
                                                                       21378, 7161, 21504, 1380,
                                                                       1408, 8442, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23436, 0, 3,
                                                                       21504, 7224, 21630, 1408,
                                                                       1436, 8526, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23604, 0, 3,
                                                                       21630, 7287, 21756, 1436,
                                                                       1464, 8610, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23772, 0, 3,
                                                                       21756, 7350, 21882, 1464,
                                                                       1492, 8694, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23940, 0, 3,
                                                                       21882, 7413, 22008, 1492,
                                                                       1520, 8778, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24108, 0, 3,
                                                                       22008, 7476, 22134, 1520,
                                                                       1548, 8862, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 24276, 0, 3,
                                                                       22260, 7770, 22428, 1604,
                                                                       1640, 9162, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 24492, 0, 3,
                                                                       22428, 7854, 22596, 1640,
                                                                       1676, 9270, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 24708, 0, 3,
                                                                       22596, 7938, 22764, 1676,
                                                                       1712, 9378, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 24924, 0, 3,
                                                                       22764, 8022, 22932, 1712,
                                                                       1748, 9486, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25140, 0, 3,
                                                                       22932, 8106, 23100, 1748,
                                                                       1784, 9594, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25356, 0, 3,
                                                                       23268, 8442, 23436, 1856,
                                                                       1892, 9918, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25572, 0, 3,
                                                                       23436, 8526, 23604, 1892,
                                                                       1928, 10026, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25788, 0, 3,
                                                                       23604, 8610, 23772, 1928,
                                                                       1964, 10134, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26004, 0, 3,
                                                                       23772, 8694, 23940, 1964,
                                                                       2000, 10242, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26220, 0, 3,
                                                                       23940, 8778, 24108, 2000,
                                                                       2036, 10350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26436, 0, 3,
                                                                       24276, 9162, 24492, 2108,
                                                                       2153, 10728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26706, 0, 3,
                                                                       24492, 9270, 24708, 2153,
                                                                       2198, 10863, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26976, 0, 3,
                                                                       24708, 9378, 24924, 2198,
                                                                       2243, 10998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 27246, 0, 3,
                                                                       24924, 9486, 25140, 2243,
                                                                       2288, 11133, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 27516, 0, 3,
                                                                       25356, 9918, 25572, 2378,
                                                                       2423, 11538, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 27786, 0, 3,
                                                                       25572, 10026, 25788, 2423,
                                                                       2468, 11673, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 28056, 0, 3,
                                                                       25788, 10134, 26004, 2468,
                                                                       2513, 11808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 28326, 0, 3,
                                                                       26004, 10242, 26220, 2513,
                                                                       2558, 11943, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 28596, 0, 3,
                                                                       26436, 10728, 26706, 2648,
                                                                       2703, 12408, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 28926, 0, 3,
                                                                       26706, 10863, 26976, 2703,
                                                                       2758, 12573, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 29256, 0, 3,
                                                                       26976, 10998, 27246, 2758,
                                                                       2813, 12738, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 29586, 0, 3,
                                                                       27516, 11538, 27786, 2923,
                                                                       2978, 13233, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 29916, 0, 3,
                                                                       27786, 11673, 28056, 2978,
                                                                       3033, 13398, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 30246, 0, 3,
                                                                       28056, 11808, 28326, 3033,
                                                                       3088, 13563, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 30576, 0, 3,
                                                                       28596, 12408, 28926, 3198,
                                                                       3264, 14124, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 30972, 0, 3,
                                                                       28926, 12573, 29256, 3264,
                                                                       3330, 14322, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 31368, 0, 3,
                                                                       29586, 13233, 29916, 3462,
                                                                       3528, 14916, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 31764, 0, 3,
                                                                       29916, 13398, 30246, 3528,
                                                                       3594, 15114, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 32160, 0, 3,
                                                                       30576, 14124, 30972, 3726,
                                                                       3804, 15780, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 32628, 0, 3,
                                                                       31368, 14916, 31764, 3960,
                                                                       4038, 16482, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33096, 3, 4194,
                                                                       4197, 16716, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33106, 3, 4197,
                                                                       4200, 16722, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33116, 3, 4200,
                                                                       4203, 16728, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33126, 3, 4203,
                                                                       4206, 16734, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33136, 3, 4206,
                                                                       4209, 16740, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33146, 3, 4209,
                                                                       4212, 16746, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33156, 3, 4212,
                                                                       4215, 16752, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33166, 3, 4215,
                                                                       4218, 16758, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33176, 3, 4218,
                                                                       4221, 16764, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33186, 3, 4221,
                                                                       4224, 16770, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33196, 3, 4224,
                                                                       4227, 16776, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33206, 3, 4227,
                                                                       4230, 16782, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33216, 3, 4236,
                                                                       4239, 16788, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33226, 3, 4239,
                                                                       4242, 16794, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33236, 3, 4242,
                                                                       4245, 16800, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33246, 3, 4245,
                                                                       4248, 16806, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33256, 3, 4248,
                                                                       4251, 16812, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33266, 3, 4251,
                                                                       4254, 16818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33276, 3, 4254,
                                                                       4257, 16824, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33286, 3, 4257,
                                                                       4260, 16830, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33296, 3, 4260,
                                                                       4263, 16836, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33306, 3, 4263,
                                                                       4266, 16842, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33316, 3, 4266,
                                                                       4269, 16848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 33326, 3, 4269,
                                                                       4272, 16854, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33336, 0, 3,
                                                                       33096, 16716, 33106,
                                                                       16860, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33366, 0, 3,
                                                                       33106, 16722, 33116,
                                                                       16878, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33396, 0, 3,
                                                                       33116, 16728, 33126,
                                                                       16896, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33426, 0, 3,
                                                                       33126, 16734, 33136,
                                                                       16914, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33456, 0, 3,
                                                                       33136, 16740, 33146,
                                                                       16932, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33486, 0, 3,
                                                                       33146, 16746, 33156,
                                                                       16950, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33516, 0, 3,
                                                                       33156, 16752, 33166,
                                                                       16968, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33546, 0, 3,
                                                                       33166, 16758, 33176,
                                                                       16986, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33576, 0, 3,
                                                                       33176, 16764, 33186,
                                                                       17004, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33606, 0, 3,
                                                                       33186, 16770, 33196,
                                                                       17022, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33636, 0, 3,
                                                                       33196, 16776, 33206,
                                                                       17040, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33666, 0, 3,
                                                                       33216, 16788, 33226,
                                                                       17058, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33696, 0, 3,
                                                                       33226, 16794, 33236,
                                                                       17076, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33726, 0, 3,
                                                                       33236, 16800, 33246,
                                                                       17094, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33756, 0, 3,
                                                                       33246, 16806, 33256,
                                                                       17112, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33786, 0, 3,
                                                                       33256, 16812, 33266,
                                                                       17130, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33816, 0, 3,
                                                                       33266, 16818, 33276,
                                                                       17148, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33846, 0, 3,
                                                                       33276, 16824, 33286,
                                                                       17166, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33876, 0, 3,
                                                                       33286, 16830, 33296,
                                                                       17184, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33906, 0, 3,
                                                                       33296, 16836, 33306,
                                                                       17202, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33936, 0, 3,
                                                                       33306, 16842, 33316,
                                                                       17220, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 33966, 0, 3,
                                                                       33316, 16848, 33326,
                                                                       17238, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 33996, 0, 3,
                                                                       33336, 16860, 33366, 4476,
                                                                       4494, 17256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34056, 0, 3,
                                                                       33366, 16878, 33396, 4494,
                                                                       4512, 17292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34116, 0, 3,
                                                                       33396, 16896, 33426, 4512,
                                                                       4530, 17328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34176, 0, 3,
                                                                       33426, 16914, 33456, 4530,
                                                                       4548, 17364, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34236, 0, 3,
                                                                       33456, 16932, 33486, 4548,
                                                                       4566, 17400, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34296, 0, 3,
                                                                       33486, 16950, 33516, 4566,
                                                                       4584, 17436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34356, 0, 3,
                                                                       33516, 16968, 33546, 4584,
                                                                       4602, 17472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34416, 0, 3,
                                                                       33546, 16986, 33576, 4602,
                                                                       4620, 17508, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34476, 0, 3,
                                                                       33576, 17004, 33606, 4620,
                                                                       4638, 17544, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34536, 0, 3,
                                                                       33606, 17022, 33636, 4638,
                                                                       4656, 17580, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34596, 0, 3,
                                                                       33666, 17058, 33696, 4692,
                                                                       4710, 17616, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34656, 0, 3,
                                                                       33696, 17076, 33726, 4710,
                                                                       4728, 17652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34716, 0, 3,
                                                                       33726, 17094, 33756, 4728,
                                                                       4746, 17688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34776, 0, 3,
                                                                       33756, 17112, 33786, 4746,
                                                                       4764, 17724, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34836, 0, 3,
                                                                       33786, 17130, 33816, 4764,
                                                                       4782, 17760, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34896, 0, 3,
                                                                       33816, 17148, 33846, 4782,
                                                                       4800, 17796, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 34956, 0, 3,
                                                                       33846, 17166, 33876, 4800,
                                                                       4818, 17832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 35016, 0, 3,
                                                                       33876, 17184, 33906, 4818,
                                                                       4836, 17868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 35076, 0, 3,
                                                                       33906, 17202, 33936, 4836,
                                                                       4854, 17904, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 35136, 0, 3,
                                                                       33936, 17220, 33966, 4854,
                                                                       4872, 17940, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35196, 0, 3,
                                                                       33996, 17256, 34056, 4908,
                                                                       4938, 17976, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35296, 0, 3,
                                                                       34056, 17292, 34116, 4938,
                                                                       4968, 18036, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35396, 0, 3,
                                                                       34116, 17328, 34176, 4968,
                                                                       4998, 18096, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35496, 0, 3,
                                                                       34176, 17364, 34236, 4998,
                                                                       5028, 18156, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35596, 0, 3,
                                                                       34236, 17400, 34296, 5028,
                                                                       5058, 18216, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35696, 0, 3,
                                                                       34296, 17436, 34356, 5058,
                                                                       5088, 18276, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35796, 0, 3,
                                                                       34356, 17472, 34416, 5088,
                                                                       5118, 18336, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35896, 0, 3,
                                                                       34416, 17508, 34476, 5118,
                                                                       5148, 18396, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 35996, 0, 3,
                                                                       34476, 17544, 34536, 5148,
                                                                       5178, 18456, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36096, 0, 3,
                                                                       34596, 17616, 34656, 5238,
                                                                       5268, 18516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36196, 0, 3,
                                                                       34656, 17652, 34716, 5268,
                                                                       5298, 18576, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36296, 0, 3,
                                                                       34716, 17688, 34776, 5298,
                                                                       5328, 18636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36396, 0, 3,
                                                                       34776, 17724, 34836, 5328,
                                                                       5358, 18696, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36496, 0, 3,
                                                                       34836, 17760, 34896, 5358,
                                                                       5388, 18756, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36596, 0, 3,
                                                                       34896, 17796, 34956, 5388,
                                                                       5418, 18816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36696, 0, 3,
                                                                       34956, 17832, 35016, 5418,
                                                                       5448, 18876, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36796, 0, 3,
                                                                       35016, 17868, 35076, 5448,
                                                                       5478, 18936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 36896, 0, 3,
                                                                       35076, 17904, 35136, 5478,
                                                                       5508, 18996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 36996, 0, 3,
                                                                       35196, 17976, 35296, 5568,
                                                                       5613, 19056, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37146, 0, 3,
                                                                       35296, 18036, 35396, 5613,
                                                                       5658, 19146, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37296, 0, 3,
                                                                       35396, 18096, 35496, 5658,
                                                                       5703, 19236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37446, 0, 3,
                                                                       35496, 18156, 35596, 5703,
                                                                       5748, 19326, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37596, 0, 3,
                                                                       35596, 18216, 35696, 5748,
                                                                       5793, 19416, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37746, 0, 3,
                                                                       35696, 18276, 35796, 5793,
                                                                       5838, 19506, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 37896, 0, 3,
                                                                       35796, 18336, 35896, 5838,
                                                                       5883, 19596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38046, 0, 3,
                                                                       35896, 18396, 35996, 5883,
                                                                       5928, 19686, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38196, 0, 3,
                                                                       36096, 18516, 36196, 6018,
                                                                       6063, 19776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38346, 0, 3,
                                                                       36196, 18576, 36296, 6063,
                                                                       6108, 19866, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38496, 0, 3,
                                                                       36296, 18636, 36396, 6108,
                                                                       6153, 19956, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38646, 0, 3,
                                                                       36396, 18696, 36496, 6153,
                                                                       6198, 20046, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38796, 0, 3,
                                                                       36496, 18756, 36596, 6198,
                                                                       6243, 20136, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38946, 0, 3,
                                                                       36596, 18816, 36696, 6243,
                                                                       6288, 20226, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39096, 0, 3,
                                                                       36696, 18876, 36796, 6288,
                                                                       6333, 20316, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39246, 0, 3,
                                                                       36796, 18936, 36896, 6333,
                                                                       6378, 20406, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 39396, 0, 3,
                                                                       36996, 19056, 37146, 6468,
                                                                       6531, 20496, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 39606, 0, 3,
                                                                       37146, 19146, 37296, 6531,
                                                                       6594, 20622, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 39816, 0, 3,
                                                                       37296, 19236, 37446, 6594,
                                                                       6657, 20748, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40026, 0, 3,
                                                                       37446, 19326, 37596, 6657,
                                                                       6720, 20874, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40236, 0, 3,
                                                                       37596, 19416, 37746, 6720,
                                                                       6783, 21000, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40446, 0, 3,
                                                                       37746, 19506, 37896, 6783,
                                                                       6846, 21126, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40656, 0, 3,
                                                                       37896, 19596, 38046, 6846,
                                                                       6909, 21252, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40866, 0, 3,
                                                                       38196, 19776, 38346, 7035,
                                                                       7098, 21378, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41076, 0, 3,
                                                                       38346, 19866, 38496, 7098,
                                                                       7161, 21504, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41286, 0, 3,
                                                                       38496, 19956, 38646, 7161,
                                                                       7224, 21630, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41496, 0, 3,
                                                                       38646, 20046, 38796, 7224,
                                                                       7287, 21756, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41706, 0, 3,
                                                                       38796, 20136, 38946, 7287,
                                                                       7350, 21882, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41916, 0, 3,
                                                                       38946, 20226, 39096, 7350,
                                                                       7413, 22008, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42126, 0, 3,
                                                                       39096, 20316, 39246, 7413,
                                                                       7476, 22134, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 42336, 0, 3,
                                                                       39396, 20496, 39606, 7602,
                                                                       7686, 22260, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 42616, 0, 3,
                                                                       39606, 20622, 39816, 7686,
                                                                       7770, 22428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 42896, 0, 3,
                                                                       39816, 20748, 40026, 7770,
                                                                       7854, 22596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43176, 0, 3,
                                                                       40026, 20874, 40236, 7854,
                                                                       7938, 22764, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43456, 0, 3,
                                                                       40236, 21000, 40446, 7938,
                                                                       8022, 22932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43736, 0, 3,
                                                                       40446, 21126, 40656, 8022,
                                                                       8106, 23100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44016, 0, 3,
                                                                       40866, 21378, 41076, 8274,
                                                                       8358, 23268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44296, 0, 3,
                                                                       41076, 21504, 41286, 8358,
                                                                       8442, 23436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44576, 0, 3,
                                                                       41286, 21630, 41496, 8442,
                                                                       8526, 23604, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44856, 0, 3,
                                                                       41496, 21756, 41706, 8526,
                                                                       8610, 23772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 45136, 0, 3,
                                                                       41706, 21882, 41916, 8610,
                                                                       8694, 23940, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 45416, 0, 3,
                                                                       41916, 22008, 42126, 8694,
                                                                       8778, 24108, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 45696, 0, 3,
                                                                       42336, 22260, 42616, 8946,
                                                                       9054, 24276, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46056, 0, 3,
                                                                       42616, 22428, 42896, 9054,
                                                                       9162, 24492, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46416, 0, 3,
                                                                       42896, 22596, 43176, 9162,
                                                                       9270, 24708, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46776, 0, 3,
                                                                       43176, 22764, 43456, 9270,
                                                                       9378, 24924, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 47136, 0, 3,
                                                                       43456, 22932, 43736, 9378,
                                                                       9486, 25140, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 47496, 0, 3,
                                                                       44016, 23268, 44296, 9702,
                                                                       9810, 25356, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 47856, 0, 3,
                                                                       44296, 23436, 44576, 9810,
                                                                       9918, 25572, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48216, 0, 3,
                                                                       44576, 23604, 44856, 9918,
                                                                       10026, 25788, ncols,
                                                                       gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48576, 0, 3,
                                                                       44856, 23772, 45136,
                                                                       10026, 10134, 26004,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48936, 0, 3,
                                                                       45136, 23940, 45416,
                                                                       10134, 10242, 26220,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 49296, 0, 3,
                                                                       45696, 24276, 46056,
                                                                       10458, 10593, 26436,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 49746, 0, 3,
                                                                       46056, 24492, 46416,
                                                                       10593, 10728, 26706,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 50196, 0, 3,
                                                                       46416, 24708, 46776,
                                                                       10728, 10863, 26976,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 50646, 0, 3,
                                                                       46776, 24924, 47136,
                                                                       10863, 10998, 27246,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 51096, 0, 3,
                                                                       47496, 25356, 47856,
                                                                       11268, 11403, 27516,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 51546, 0, 3,
                                                                       47856, 25572, 48216,
                                                                       11403, 11538, 27786,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 51996, 0, 3,
                                                                       48216, 25788, 48576,
                                                                       11538, 11673, 28056,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 52446, 0, 3,
                                                                       48576, 26004, 48936,
                                                                       11673, 11808, 28326,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 52896, 0, 3,
                                                                       49296, 26436, 49746,
                                                                       12078, 12243, 28596,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 53446, 0, 3,
                                                                       49746, 26706, 50196,
                                                                       12243, 12408, 28926,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 53996, 0, 3,
                                                                       50196, 26976, 50646,
                                                                       12408, 12573, 29256,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 54546, 0, 3,
                                                                       51096, 27516, 51546,
                                                                       12903, 13068, 29586,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 55096, 0, 3,
                                                                       51546, 27786, 51996,
                                                                       13068, 13233, 29916,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 55646, 0, 3,
                                                                       51996, 28056, 52446,
                                                                       13233, 13398, 30246,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 56196, 0, 3,
                                                                       52896, 28596, 53446,
                                                                       13728, 13926, 30576,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 56856, 0, 3,
                                                                       53446, 28926, 53996,
                                                                       13926, 14124, 30972,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 57516, 0, 3,
                                                                       54546, 29586, 55096,
                                                                       14520, 14718, 31368,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 58176, 0, 3,
                                                                       55096, 29916, 55646,
                                                                       14718, 14916, 31764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 58836, 0, 3,
                                                                       56196, 30576, 56856,
                                                                       15312, 15546, 32160,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 59616, 0, 3,
                                                                       57516, 31368, 58176,
                                                                       16014, 16248, 32628,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 60396, 42336, 280, ncols);

                    simdfunc::contract_primitives(buffer, 60872, 44016, 280, ncols);

                    simdfunc::contract_primitives(buffer, 61348, 45696, 360, ncols);

                    simdfunc::contract_primitives(buffer, 61960, 47496, 360, ncols);

                    simdfunc::contract_primitives(buffer, 62572, 49296, 450, ncols);

                    simdfunc::contract_primitives(buffer, 63337, 51096, 450, ncols);

                    simdfunc::contract_primitives(buffer, 64102, 52896, 550, ncols);

                    simdfunc::contract_primitives(buffer, 65037, 54546, 550, ncols);

                    simdfunc::contract_primitives(buffer, 65972, 56196, 660, ncols);

                    simdfunc::contract_primitives(buffer, 67094, 57516, 660, ncols);

                    simdfunc::contract_primitives(buffer, 68216, 58836, 780, ncols);

                    simdfunc::contract_primitives(buffer, 69542, 59616, 780, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 60676, 60396, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 61152, 60872, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 61708, 61348, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 62320, 61960, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 63022, 62572, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 63787, 63337, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 64652, 64102, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 65587, 65037, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 66632, 65972, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 67754, 67094, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 68996, 68216, 78, 1, nmax);

        simdtrf::transform_f_inner(buffer, 70322, 69542, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 70868, 60676, 61708, 7, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 71456, 61152, 62320, 7, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 72044, 61708, 63022, 7, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 72800, 62320, 63787, 7, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 73556, 63022, 64652, 7, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 74501, 63787, 65587, 7, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 75446, 64652, 66632, 7, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 76601, 65587, 67754, 7, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 77756, 66632, 68996, 7, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 79142, 67754, 70322, 7, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 80528, 70868, 72044, 7, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 81704, 71456, 72800, 7, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 82880, 72044, 73556, 7, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 84392, 72800, 74501, 7, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 85904, 73556, 75446, 7, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 87794, 74501, 76601, 7, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 89684, 75446, 77756, 7, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 91994, 76601, 79142, 7, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 94304, 80528, 82880, 7, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 96264, 81704, 84392, 7, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 98224, 82880, 85904, 7, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 100744, 84392, 87794, 7, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 103264, 85904, 89684, 7, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 106414, 87794, 91994, 7, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 109564, 94304, 98224, 7, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 112504, 96264, 100744, 7, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 115444, 98224, 103264, 7, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 119224, 100744, 106414, 7, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 123004, 109564, 115444, 7, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 127120, 112504, 119224, 7, nmax);

        simdtrf::transform_i_inner(buffer, 131236, 127120, 21, 7, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 131236, 91, nmax);

        simdtrf::transform_i_inner(buffer, 131236, 123004, 21, 7, nmax);

        simdtrf::transform_h_outer(values + 1001 * nvalues + n * npairs, nvalues, buffer, 131236,
                                   91, nmax);
    }

    for (size_t m = 0; m < 2002; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
