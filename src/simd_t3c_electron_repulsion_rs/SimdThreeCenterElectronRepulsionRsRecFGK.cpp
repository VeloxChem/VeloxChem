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


#include "SimdThreeCenterElectronRepulsionRsRecFGK.hpp"

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
#include "SimdTransferDG.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferFG.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_fgk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_fgk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 145862, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1890 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 145862, 117572, 9660, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2108, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2111, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2114, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2117, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2120, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2123, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2126, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2129, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2132, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2135, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2138, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2141, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2144, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2147, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2150, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2153, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2156, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2159, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2162, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2165, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2168, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2171, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2174, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2177, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2180, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2183, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2186, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2189, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2192, 3, 9, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2201, 3, 10, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2210, 3, 11, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2219, 3, 12, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2228, 3, 13, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2237, 3, 14, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2246, 3, 15, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2255, 3, 16, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2264, 3, 17, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2273, 3, 18, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2282, 3, 19, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2291, 3, 24, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2300, 3, 25, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2309, 3, 26, 87,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2318, 3, 27, 90,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2327, 3, 28, 93,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2336, 3, 29, 96,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2345, 3, 30, 99,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2354, 3, 31, 102,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2363, 3, 32, 105,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2372, 3, 33, 108,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2381, 3, 34, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2390, 3, 36, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2408, 3, 39, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2426, 3, 42, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2444, 3, 45, 132,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2462, 3, 48, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2480, 3, 51, 144,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2498, 3, 54, 150,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2516, 3, 57, 156,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2534, 3, 60, 162,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2552, 3, 63, 168,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2570, 3, 66, 174,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2588, 3, 69, 180,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2606, 3, 75, 186,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2624, 3, 78, 192,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2642, 3, 81, 198,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2660, 3, 84, 204,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2678, 3, 87, 210,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2696, 3, 90, 216,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2714, 3, 93, 222,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2732, 3, 96, 228,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2750, 3, 99, 234,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2768, 3, 102, 240,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2786, 3, 105, 246,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2804, 3, 108, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2822, 3, 114, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2852, 3, 120, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2882, 3, 126, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2912, 3, 132, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2942, 3, 138, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2972, 3, 144, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3002, 3, 150, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3032, 3, 156, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3062, 3, 162, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3092, 3, 168, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3122, 3, 174, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3152, 3, 186, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3182, 3, 192, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3212, 3, 198, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3242, 3, 204, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3272, 3, 210, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3302, 3, 216, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3332, 3, 222, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3362, 3, 228, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3392, 3, 234, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3422, 3, 240, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3452, 3, 246, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3482, 3, 258, 478,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3527, 3, 268, 493,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3572, 3, 278, 508,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3617, 3, 288, 523,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3662, 3, 298, 538,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3707, 3, 308, 553,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3752, 3, 318, 568,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3797, 3, 328, 583,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3842, 3, 338, 598,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3887, 3, 348, 613,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3932, 3, 368, 628,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3977, 3, 378, 643,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4022, 3, 388, 658,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4067, 3, 398, 673,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4112, 3, 408, 688,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4157, 3, 418, 703,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4202, 3, 428, 718,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4247, 3, 438, 733,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4292, 3, 448, 748,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4337, 3, 458, 763,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4382, 3, 478, 778,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4445, 3, 493, 799,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4508, 3, 508, 820,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4571, 3, 523, 841,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4634, 3, 538, 862,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4697, 3, 553, 883,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4760, 3, 568, 904,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4823, 3, 583, 925,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4886, 3, 598, 946,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4949, 3, 628, 967,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5012, 3, 643, 988,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5075, 3, 658,
                                                                       1009, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5138, 3, 673,
                                                                       1030, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5201, 3, 688,
                                                                       1051, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5264, 3, 703,
                                                                       1072, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5327, 3, 718,
                                                                       1093, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5390, 3, 733,
                                                                       1114, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5453, 3, 748,
                                                                       1135, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5516, 3, 778,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5600, 3, 799,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5684, 3, 820,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5768, 3, 841,
                                                                       1240, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5852, 3, 862,
                                                                       1268, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5936, 3, 883,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6020, 3, 904,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6104, 3, 925,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6188, 3, 967,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6272, 3, 988,
                                                                       1408, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6356, 3, 1009,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6440, 3, 1030,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6524, 3, 1051,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6608, 3, 1072,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6692, 3, 1093,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6776, 3, 1114,
                                                                       1576, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6860, 3, 1156,
                                                                       1604, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6968, 3, 1184,
                                                                       1640, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7076, 3, 1212,
                                                                       1676, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7184, 3, 1240,
                                                                       1712, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7292, 3, 1268,
                                                                       1748, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7400, 3, 1296,
                                                                       1784, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7508, 3, 1324,
                                                                       1820, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7616, 3, 1380,
                                                                       1856, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7724, 3, 1408,
                                                                       1892, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7832, 3, 1436,
                                                                       1928, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7940, 3, 1464,
                                                                       1964, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8048, 3, 1492,
                                                                       2000, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8156, 3, 1520,
                                                                       2036, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8264, 3, 1548,
                                                                       2072, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8372, 3, 7, 8,
                                                                       2114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8378, 3, 8, 9,
                                                                       2117, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8384, 3, 9, 10,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8390, 3, 10, 11,
                                                                       2123, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8396, 3, 11, 12,
                                                                       2126, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8402, 3, 12, 13,
                                                                       2129, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8408, 3, 13, 14,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8414, 3, 14, 15,
                                                                       2135, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8420, 3, 15, 16,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8426, 3, 16, 17,
                                                                       2141, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8432, 3, 17, 18,
                                                                       2144, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8438, 3, 18, 19,
                                                                       2147, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8444, 3, 22, 23,
                                                                       2156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8450, 3, 23, 24,
                                                                       2159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8456, 3, 24, 25,
                                                                       2162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8462, 3, 25, 26,
                                                                       2165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8468, 3, 26, 27,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8474, 3, 27, 28,
                                                                       2171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8480, 3, 28, 29,
                                                                       2174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8486, 3, 29, 30,
                                                                       2177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8492, 3, 30, 31,
                                                                       2180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8498, 3, 31, 32,
                                                                       2183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8504, 3, 32, 33,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8510, 3, 33, 34,
                                                                       2189, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8516, 0, 3, 8372,
                                                                       2114, 8378, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8534, 0, 3, 8378,
                                                                       2117, 8384, 2201, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8552, 0, 3, 8384,
                                                                       2120, 8390, 2210, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8570, 0, 3, 8390,
                                                                       2123, 8396, 2219, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8588, 0, 3, 8396,
                                                                       2126, 8402, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8606, 0, 3, 8402,
                                                                       2129, 8408, 2237, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8624, 0, 3, 8408,
                                                                       2132, 8414, 2246, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8642, 0, 3, 8414,
                                                                       2135, 8420, 2255, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8660, 0, 3, 8420,
                                                                       2138, 8426, 2264, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8678, 0, 3, 8426,
                                                                       2141, 8432, 2273, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8696, 0, 3, 8432,
                                                                       2144, 8438, 2282, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8714, 0, 3, 8444,
                                                                       2156, 8450, 2291, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8732, 0, 3, 8450,
                                                                       2159, 8456, 2300, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8750, 0, 3, 8456,
                                                                       2162, 8462, 2309, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8768, 0, 3, 8462,
                                                                       2165, 8468, 2318, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8786, 0, 3, 8468,
                                                                       2168, 8474, 2327, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8804, 0, 3, 8474,
                                                                       2171, 8480, 2336, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8822, 0, 3, 8480,
                                                                       2174, 8486, 2345, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8840, 0, 3, 8486,
                                                                       2177, 8492, 2354, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8858, 0, 3, 8492,
                                                                       2180, 8498, 2363, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8876, 0, 3, 8498,
                                                                       2183, 8504, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8894, 0, 3, 8504,
                                                                       2186, 8510, 2381, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8912, 0, 3, 8516,
                                                                       2192, 8534, 114, 120,
                                                                       2426, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8948, 0, 3, 8534,
                                                                       2201, 8552, 120, 126,
                                                                       2444, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8984, 0, 3, 8552,
                                                                       2210, 8570, 126, 132,
                                                                       2462, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9020, 0, 3, 8570,
                                                                       2219, 8588, 132, 138,
                                                                       2480, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9056, 0, 3, 8588,
                                                                       2228, 8606, 138, 144,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9092, 0, 3, 8606,
                                                                       2237, 8624, 144, 150,
                                                                       2516, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9128, 0, 3, 8624,
                                                                       2246, 8642, 150, 156,
                                                                       2534, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9164, 0, 3, 8642,
                                                                       2255, 8660, 156, 162,
                                                                       2552, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9200, 0, 3, 8660,
                                                                       2264, 8678, 162, 168,
                                                                       2570, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9236, 0, 3, 8678,
                                                                       2273, 8696, 168, 174,
                                                                       2588, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9272, 0, 3, 8714,
                                                                       2291, 8732, 186, 192,
                                                                       2642, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9308, 0, 3, 8732,
                                                                       2300, 8750, 192, 198,
                                                                       2660, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9344, 0, 3, 8750,
                                                                       2309, 8768, 198, 204,
                                                                       2678, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9380, 0, 3, 8768,
                                                                       2318, 8786, 204, 210,
                                                                       2696, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9416, 0, 3, 8786,
                                                                       2327, 8804, 210, 216,
                                                                       2714, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9452, 0, 3, 8804,
                                                                       2336, 8822, 216, 222,
                                                                       2732, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9488, 0, 3, 8822,
                                                                       2345, 8840, 222, 228,
                                                                       2750, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9524, 0, 3, 8840,
                                                                       2354, 8858, 228, 234,
                                                                       2768, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9560, 0, 3, 8858,
                                                                       2363, 8876, 234, 240,
                                                                       2786, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9596, 0, 3, 8876,
                                                                       2372, 8894, 240, 246,
                                                                       2804, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9632, 0, 3, 8912,
                                                                       2426, 8948, 258, 268,
                                                                       2882, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9692, 0, 3, 8948,
                                                                       2444, 8984, 268, 278,
                                                                       2912, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9752, 0, 3, 8984,
                                                                       2462, 9020, 278, 288,
                                                                       2942, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9812, 0, 3, 9020,
                                                                       2480, 9056, 288, 298,
                                                                       2972, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9872, 0, 3, 9056,
                                                                       2498, 9092, 298, 308,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9932, 0, 3, 9092,
                                                                       2516, 9128, 308, 318,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9992, 0, 3, 9128,
                                                                       2534, 9164, 318, 328,
                                                                       3062, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10052, 0, 3, 9164,
                                                                       2552, 9200, 328, 338,
                                                                       3092, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10112, 0, 3, 9200,
                                                                       2570, 9236, 338, 348,
                                                                       3122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10172, 0, 3, 9272,
                                                                       2642, 9308, 368, 378,
                                                                       3212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10232, 0, 3, 9308,
                                                                       2660, 9344, 378, 388,
                                                                       3242, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10292, 0, 3, 9344,
                                                                       2678, 9380, 388, 398,
                                                                       3272, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10352, 0, 3, 9380,
                                                                       2696, 9416, 398, 408,
                                                                       3302, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10412, 0, 3, 9416,
                                                                       2714, 9452, 408, 418,
                                                                       3332, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10472, 0, 3, 9452,
                                                                       2732, 9488, 418, 428,
                                                                       3362, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10532, 0, 3, 9488,
                                                                       2750, 9524, 428, 438,
                                                                       3392, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10592, 0, 3, 9524,
                                                                       2768, 9560, 438, 448,
                                                                       3422, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10652, 0, 3, 9560,
                                                                       2786, 9596, 448, 458,
                                                                       3452, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10712, 0, 3, 9632,
                                                                       2882, 9692, 478, 493,
                                                                       3572, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10802, 0, 3, 9692,
                                                                       2912, 9752, 493, 508,
                                                                       3617, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10892, 0, 3, 9752,
                                                                       2942, 9812, 508, 523,
                                                                       3662, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10982, 0, 3, 9812,
                                                                       2972, 9872, 523, 538,
                                                                       3707, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11072, 0, 3, 9872,
                                                                       3002, 9932, 538, 553,
                                                                       3752, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11162, 0, 3, 9932,
                                                                       3032, 9992, 553, 568,
                                                                       3797, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11252, 0, 3, 9992,
                                                                       3062, 10052, 568, 583,
                                                                       3842, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11342, 0, 3,
                                                                       10052, 3092, 10112, 583,
                                                                       598, 3887, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11432, 0, 3,
                                                                       10172, 3212, 10232, 628,
                                                                       643, 4022, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11522, 0, 3,
                                                                       10232, 3242, 10292, 643,
                                                                       658, 4067, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11612, 0, 3,
                                                                       10292, 3272, 10352, 658,
                                                                       673, 4112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11702, 0, 3,
                                                                       10352, 3302, 10412, 673,
                                                                       688, 4157, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11792, 0, 3,
                                                                       10412, 3332, 10472, 688,
                                                                       703, 4202, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11882, 0, 3,
                                                                       10472, 3362, 10532, 703,
                                                                       718, 4247, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11972, 0, 3,
                                                                       10532, 3392, 10592, 718,
                                                                       733, 4292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12062, 0, 3,
                                                                       10592, 3422, 10652, 733,
                                                                       748, 4337, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12152, 0, 3,
                                                                       10712, 3572, 10802, 778,
                                                                       799, 4508, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12278, 0, 3,
                                                                       10802, 3617, 10892, 799,
                                                                       820, 4571, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12404, 0, 3,
                                                                       10892, 3662, 10982, 820,
                                                                       841, 4634, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12530, 0, 3,
                                                                       10982, 3707, 11072, 841,
                                                                       862, 4697, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12656, 0, 3,
                                                                       11072, 3752, 11162, 862,
                                                                       883, 4760, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12782, 0, 3,
                                                                       11162, 3797, 11252, 883,
                                                                       904, 4823, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12908, 0, 3,
                                                                       11252, 3842, 11342, 904,
                                                                       925, 4886, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13034, 0, 3,
                                                                       11432, 4022, 11522, 967,
                                                                       988, 5075, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13160, 0, 3,
                                                                       11522, 4067, 11612, 988,
                                                                       1009, 5138, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13286, 0, 3,
                                                                       11612, 4112, 11702, 1009,
                                                                       1030, 5201, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13412, 0, 3,
                                                                       11702, 4157, 11792, 1030,
                                                                       1051, 5264, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13538, 0, 3,
                                                                       11792, 4202, 11882, 1051,
                                                                       1072, 5327, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13664, 0, 3,
                                                                       11882, 4247, 11972, 1072,
                                                                       1093, 5390, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13790, 0, 3,
                                                                       11972, 4292, 12062, 1093,
                                                                       1114, 5453, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13916, 0, 3,
                                                                       12152, 4508, 12278, 1156,
                                                                       1184, 5684, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14084, 0, 3,
                                                                       12278, 4571, 12404, 1184,
                                                                       1212, 5768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14252, 0, 3,
                                                                       12404, 4634, 12530, 1212,
                                                                       1240, 5852, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14420, 0, 3,
                                                                       12530, 4697, 12656, 1240,
                                                                       1268, 5936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14588, 0, 3,
                                                                       12656, 4760, 12782, 1268,
                                                                       1296, 6020, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14756, 0, 3,
                                                                       12782, 4823, 12908, 1296,
                                                                       1324, 6104, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14924, 0, 3,
                                                                       13034, 5075, 13160, 1380,
                                                                       1408, 6356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15092, 0, 3,
                                                                       13160, 5138, 13286, 1408,
                                                                       1436, 6440, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15260, 0, 3,
                                                                       13286, 5201, 13412, 1436,
                                                                       1464, 6524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15428, 0, 3,
                                                                       13412, 5264, 13538, 1464,
                                                                       1492, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15596, 0, 3,
                                                                       13538, 5327, 13664, 1492,
                                                                       1520, 6692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15764, 0, 3,
                                                                       13664, 5390, 13790, 1520,
                                                                       1548, 6776, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15932, 0, 3,
                                                                       13916, 5684, 14084, 1604,
                                                                       1640, 7076, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16148, 0, 3,
                                                                       14084, 5768, 14252, 1640,
                                                                       1676, 7184, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16364, 0, 3,
                                                                       14252, 5852, 14420, 1676,
                                                                       1712, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16580, 0, 3,
                                                                       14420, 5936, 14588, 1712,
                                                                       1748, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16796, 0, 3,
                                                                       14588, 6020, 14756, 1748,
                                                                       1784, 7508, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17012, 0, 3,
                                                                       14924, 6356, 15092, 1856,
                                                                       1892, 7832, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17228, 0, 3,
                                                                       15092, 6440, 15260, 1892,
                                                                       1928, 7940, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17444, 0, 3,
                                                                       15260, 6524, 15428, 1928,
                                                                       1964, 8048, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17660, 0, 3,
                                                                       15428, 6608, 15596, 1964,
                                                                       2000, 8156, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17876, 0, 3,
                                                                       15596, 6692, 15764, 2000,
                                                                       2036, 8264, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18092, 3, 2108,
                                                                       2111, 8372, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18102, 3, 2111,
                                                                       2114, 8378, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18112, 3, 2114,
                                                                       2117, 8384, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18122, 3, 2117,
                                                                       2120, 8390, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18132, 3, 2120,
                                                                       2123, 8396, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18142, 3, 2123,
                                                                       2126, 8402, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18152, 3, 2126,
                                                                       2129, 8408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18162, 3, 2129,
                                                                       2132, 8414, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18172, 3, 2132,
                                                                       2135, 8420, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18182, 3, 2135,
                                                                       2138, 8426, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18192, 3, 2138,
                                                                       2141, 8432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18202, 3, 2141,
                                                                       2144, 8438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18212, 3, 2150,
                                                                       2153, 8444, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18222, 3, 2153,
                                                                       2156, 8450, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18232, 3, 2156,
                                                                       2159, 8456, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18242, 3, 2159,
                                                                       2162, 8462, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18252, 3, 2162,
                                                                       2165, 8468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18262, 3, 2165,
                                                                       2168, 8474, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18272, 3, 2168,
                                                                       2171, 8480, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18282, 3, 2171,
                                                                       2174, 8486, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18292, 3, 2174,
                                                                       2177, 8492, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18302, 3, 2177,
                                                                       2180, 8498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18312, 3, 2180,
                                                                       2183, 8504, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18322, 3, 2183,
                                                                       2186, 8510, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18332, 0, 3,
                                                                       18092, 8372, 18102, 8516,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18362, 0, 3,
                                                                       18102, 8378, 18112, 8534,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18392, 0, 3,
                                                                       18112, 8384, 18122, 8552,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18422, 0, 3,
                                                                       18122, 8390, 18132, 8570,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18452, 0, 3,
                                                                       18132, 8396, 18142, 8588,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18482, 0, 3,
                                                                       18142, 8402, 18152, 8606,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18512, 0, 3,
                                                                       18152, 8408, 18162, 8624,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18542, 0, 3,
                                                                       18162, 8414, 18172, 8642,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18572, 0, 3,
                                                                       18172, 8420, 18182, 8660,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18602, 0, 3,
                                                                       18182, 8426, 18192, 8678,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18632, 0, 3,
                                                                       18192, 8432, 18202, 8696,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18662, 0, 3,
                                                                       18212, 8444, 18222, 8714,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18692, 0, 3,
                                                                       18222, 8450, 18232, 8732,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18722, 0, 3,
                                                                       18232, 8456, 18242, 8750,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18752, 0, 3,
                                                                       18242, 8462, 18252, 8768,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18782, 0, 3,
                                                                       18252, 8468, 18262, 8786,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18812, 0, 3,
                                                                       18262, 8474, 18272, 8804,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18842, 0, 3,
                                                                       18272, 8480, 18282, 8822,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18872, 0, 3,
                                                                       18282, 8486, 18292, 8840,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18902, 0, 3,
                                                                       18292, 8492, 18302, 8858,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18932, 0, 3,
                                                                       18302, 8498, 18312, 8876,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18962, 0, 3,
                                                                       18312, 8504, 18322, 8894,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18992, 0, 3,
                                                                       18332, 8516, 18362, 2390,
                                                                       2408, 8912, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19052, 0, 3,
                                                                       18362, 8534, 18392, 2408,
                                                                       2426, 8948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19112, 0, 3,
                                                                       18392, 8552, 18422, 2426,
                                                                       2444, 8984, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19172, 0, 3,
                                                                       18422, 8570, 18452, 2444,
                                                                       2462, 9020, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19232, 0, 3,
                                                                       18452, 8588, 18482, 2462,
                                                                       2480, 9056, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19292, 0, 3,
                                                                       18482, 8606, 18512, 2480,
                                                                       2498, 9092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19352, 0, 3,
                                                                       18512, 8624, 18542, 2498,
                                                                       2516, 9128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19412, 0, 3,
                                                                       18542, 8642, 18572, 2516,
                                                                       2534, 9164, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19472, 0, 3,
                                                                       18572, 8660, 18602, 2534,
                                                                       2552, 9200, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19532, 0, 3,
                                                                       18602, 8678, 18632, 2552,
                                                                       2570, 9236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19592, 0, 3,
                                                                       18662, 8714, 18692, 2606,
                                                                       2624, 9272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19652, 0, 3,
                                                                       18692, 8732, 18722, 2624,
                                                                       2642, 9308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19712, 0, 3,
                                                                       18722, 8750, 18752, 2642,
                                                                       2660, 9344, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19772, 0, 3,
                                                                       18752, 8768, 18782, 2660,
                                                                       2678, 9380, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19832, 0, 3,
                                                                       18782, 8786, 18812, 2678,
                                                                       2696, 9416, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19892, 0, 3,
                                                                       18812, 8804, 18842, 2696,
                                                                       2714, 9452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 19952, 0, 3,
                                                                       18842, 8822, 18872, 2714,
                                                                       2732, 9488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20012, 0, 3,
                                                                       18872, 8840, 18902, 2732,
                                                                       2750, 9524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20072, 0, 3,
                                                                       18902, 8858, 18932, 2750,
                                                                       2768, 9560, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 20132, 0, 3,
                                                                       18932, 8876, 18962, 2768,
                                                                       2786, 9596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20192, 0, 3,
                                                                       18992, 8912, 19052, 2822,
                                                                       2852, 9632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20292, 0, 3,
                                                                       19052, 8948, 19112, 2852,
                                                                       2882, 9692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20392, 0, 3,
                                                                       19112, 8984, 19172, 2882,
                                                                       2912, 9752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20492, 0, 3,
                                                                       19172, 9020, 19232, 2912,
                                                                       2942, 9812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20592, 0, 3,
                                                                       19232, 9056, 19292, 2942,
                                                                       2972, 9872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20692, 0, 3,
                                                                       19292, 9092, 19352, 2972,
                                                                       3002, 9932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20792, 0, 3,
                                                                       19352, 9128, 19412, 3002,
                                                                       3032, 9992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20892, 0, 3,
                                                                       19412, 9164, 19472, 3032,
                                                                       3062, 10052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 20992, 0, 3,
                                                                       19472, 9200, 19532, 3062,
                                                                       3092, 10112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21092, 0, 3,
                                                                       19592, 9272, 19652, 3152,
                                                                       3182, 10172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21192, 0, 3,
                                                                       19652, 9308, 19712, 3182,
                                                                       3212, 10232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21292, 0, 3,
                                                                       19712, 9344, 19772, 3212,
                                                                       3242, 10292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21392, 0, 3,
                                                                       19772, 9380, 19832, 3242,
                                                                       3272, 10352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21492, 0, 3,
                                                                       19832, 9416, 19892, 3272,
                                                                       3302, 10412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21592, 0, 3,
                                                                       19892, 9452, 19952, 3302,
                                                                       3332, 10472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21692, 0, 3,
                                                                       19952, 9488, 20012, 3332,
                                                                       3362, 10532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21792, 0, 3,
                                                                       20012, 9524, 20072, 3362,
                                                                       3392, 10592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 21892, 0, 3,
                                                                       20072, 9560, 20132, 3392,
                                                                       3422, 10652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21992, 0, 3,
                                                                       20192, 9632, 20292, 3482,
                                                                       3527, 10712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 22142, 0, 3,
                                                                       20292, 9692, 20392, 3527,
                                                                       3572, 10802, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 22292, 0, 3,
                                                                       20392, 9752, 20492, 3572,
                                                                       3617, 10892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 22442, 0, 3,
                                                                       20492, 9812, 20592, 3617,
                                                                       3662, 10982, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 22592, 0, 3,
                                                                       20592, 9872, 20692, 3662,
                                                                       3707, 11072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 22742, 0, 3,
                                                                       20692, 9932, 20792, 3707,
                                                                       3752, 11162, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 22892, 0, 3,
                                                                       20792, 9992, 20892, 3752,
                                                                       3797, 11252, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23042, 0, 3,
                                                                       20892, 10052, 20992, 3797,
                                                                       3842, 11342, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23192, 0, 3,
                                                                       21092, 10172, 21192, 3932,
                                                                       3977, 11432, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23342, 0, 3,
                                                                       21192, 10232, 21292, 3977,
                                                                       4022, 11522, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23492, 0, 3,
                                                                       21292, 10292, 21392, 4022,
                                                                       4067, 11612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23642, 0, 3,
                                                                       21392, 10352, 21492, 4067,
                                                                       4112, 11702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23792, 0, 3,
                                                                       21492, 10412, 21592, 4112,
                                                                       4157, 11792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23942, 0, 3,
                                                                       21592, 10472, 21692, 4157,
                                                                       4202, 11882, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24092, 0, 3,
                                                                       21692, 10532, 21792, 4202,
                                                                       4247, 11972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24242, 0, 3,
                                                                       21792, 10592, 21892, 4247,
                                                                       4292, 12062, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 24392, 0, 3,
                                                                       21992, 10712, 22142, 4382,
                                                                       4445, 12152, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 24602, 0, 3,
                                                                       22142, 10802, 22292, 4445,
                                                                       4508, 12278, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 24812, 0, 3,
                                                                       22292, 10892, 22442, 4508,
                                                                       4571, 12404, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25022, 0, 3,
                                                                       22442, 10982, 22592, 4571,
                                                                       4634, 12530, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25232, 0, 3,
                                                                       22592, 11072, 22742, 4634,
                                                                       4697, 12656, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25442, 0, 3,
                                                                       22742, 11162, 22892, 4697,
                                                                       4760, 12782, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25652, 0, 3,
                                                                       22892, 11252, 23042, 4760,
                                                                       4823, 12908, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25862, 0, 3,
                                                                       23192, 11432, 23342, 4949,
                                                                       5012, 13034, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26072, 0, 3,
                                                                       23342, 11522, 23492, 5012,
                                                                       5075, 13160, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26282, 0, 3,
                                                                       23492, 11612, 23642, 5075,
                                                                       5138, 13286, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26492, 0, 3,
                                                                       23642, 11702, 23792, 5138,
                                                                       5201, 13412, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26702, 0, 3,
                                                                       23792, 11792, 23942, 5201,
                                                                       5264, 13538, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26912, 0, 3,
                                                                       23942, 11882, 24092, 5264,
                                                                       5327, 13664, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27122, 0, 3,
                                                                       24092, 11972, 24242, 5327,
                                                                       5390, 13790, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27332, 0, 3,
                                                                       24392, 12152, 24602, 5516,
                                                                       5600, 13916, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27612, 0, 3,
                                                                       24602, 12278, 24812, 5600,
                                                                       5684, 14084, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27892, 0, 3,
                                                                       24812, 12404, 25022, 5684,
                                                                       5768, 14252, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28172, 0, 3,
                                                                       25022, 12530, 25232, 5768,
                                                                       5852, 14420, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28452, 0, 3,
                                                                       25232, 12656, 25442, 5852,
                                                                       5936, 14588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28732, 0, 3,
                                                                       25442, 12782, 25652, 5936,
                                                                       6020, 14756, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29012, 0, 3,
                                                                       25862, 13034, 26072, 6188,
                                                                       6272, 14924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29292, 0, 3,
                                                                       26072, 13160, 26282, 6272,
                                                                       6356, 15092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29572, 0, 3,
                                                                       26282, 13286, 26492, 6356,
                                                                       6440, 15260, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29852, 0, 3,
                                                                       26492, 13412, 26702, 6440,
                                                                       6524, 15428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30132, 0, 3,
                                                                       26702, 13538, 26912, 6524,
                                                                       6608, 15596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30412, 0, 3,
                                                                       26912, 13664, 27122, 6608,
                                                                       6692, 15764, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30692, 0, 3,
                                                                       27332, 13916, 27612, 6860,
                                                                       6968, 15932, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 31052, 0, 3,
                                                                       27612, 14084, 27892, 6968,
                                                                       7076, 16148, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 31412, 0, 3,
                                                                       27892, 14252, 28172, 7076,
                                                                       7184, 16364, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 31772, 0, 3,
                                                                       28172, 14420, 28452, 7184,
                                                                       7292, 16580, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32132, 0, 3,
                                                                       28452, 14588, 28732, 7292,
                                                                       7400, 16796, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32492, 0, 3,
                                                                       29012, 14924, 29292, 7616,
                                                                       7724, 17012, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32852, 0, 3,
                                                                       29292, 15092, 29572, 7724,
                                                                       7832, 17228, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33212, 0, 3,
                                                                       29572, 15260, 29852, 7832,
                                                                       7940, 17444, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33572, 0, 3,
                                                                       29852, 15428, 30132, 7940,
                                                                       8048, 17660, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33932, 0, 3,
                                                                       30132, 15596, 30412, 8048,
                                                                       8156, 17876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34292, 3, 8372,
                                                                       8378, 18112, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34307, 3, 8378,
                                                                       8384, 18122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34322, 3, 8384,
                                                                       8390, 18132, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34337, 3, 8390,
                                                                       8396, 18142, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34352, 3, 8396,
                                                                       8402, 18152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34367, 3, 8402,
                                                                       8408, 18162, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34382, 3, 8408,
                                                                       8414, 18172, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34397, 3, 8414,
                                                                       8420, 18182, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34412, 3, 8420,
                                                                       8426, 18192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34427, 3, 8426,
                                                                       8432, 18202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34442, 3, 8444,
                                                                       8450, 18232, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34457, 3, 8450,
                                                                       8456, 18242, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34472, 3, 8456,
                                                                       8462, 18252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34487, 3, 8462,
                                                                       8468, 18262, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34502, 3, 8468,
                                                                       8474, 18272, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34517, 3, 8474,
                                                                       8480, 18282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34532, 3, 8480,
                                                                       8486, 18292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34547, 3, 8486,
                                                                       8492, 18302, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34562, 3, 8492,
                                                                       8498, 18312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 34577, 3, 8498,
                                                                       8504, 18322, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34592, 0, 3,
                                                                       34292, 18112, 34307, 8516,
                                                                       8534, 18392, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34637, 0, 3,
                                                                       34307, 18122, 34322, 8534,
                                                                       8552, 18422, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34682, 0, 3,
                                                                       34322, 18132, 34337, 8552,
                                                                       8570, 18452, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34727, 0, 3,
                                                                       34337, 18142, 34352, 8570,
                                                                       8588, 18482, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34772, 0, 3,
                                                                       34352, 18152, 34367, 8588,
                                                                       8606, 18512, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34817, 0, 3,
                                                                       34367, 18162, 34382, 8606,
                                                                       8624, 18542, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34862, 0, 3,
                                                                       34382, 18172, 34397, 8624,
                                                                       8642, 18572, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34907, 0, 3,
                                                                       34397, 18182, 34412, 8642,
                                                                       8660, 18602, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34952, 0, 3,
                                                                       34412, 18192, 34427, 8660,
                                                                       8678, 18632, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34997, 0, 3,
                                                                       34442, 18232, 34457, 8714,
                                                                       8732, 18722, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35042, 0, 3,
                                                                       34457, 18242, 34472, 8732,
                                                                       8750, 18752, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35087, 0, 3,
                                                                       34472, 18252, 34487, 8750,
                                                                       8768, 18782, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35132, 0, 3,
                                                                       34487, 18262, 34502, 8768,
                                                                       8786, 18812, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35177, 0, 3,
                                                                       34502, 18272, 34517, 8786,
                                                                       8804, 18842, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35222, 0, 3,
                                                                       34517, 18282, 34532, 8804,
                                                                       8822, 18872, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35267, 0, 3,
                                                                       34532, 18292, 34547, 8822,
                                                                       8840, 18902, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35312, 0, 3,
                                                                       34547, 18302, 34562, 8840,
                                                                       8858, 18932, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35357, 0, 3,
                                                                       34562, 18312, 34577, 8858,
                                                                       8876, 18962, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35402, 0, 3,
                                                                       34592, 18392, 34637, 8912,
                                                                       8948, 19112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35492, 0, 3,
                                                                       34637, 18422, 34682, 8948,
                                                                       8984, 19172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35582, 0, 3,
                                                                       34682, 18452, 34727, 8984,
                                                                       9020, 19232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35672, 0, 3,
                                                                       34727, 18482, 34772, 9020,
                                                                       9056, 19292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35762, 0, 3,
                                                                       34772, 18512, 34817, 9056,
                                                                       9092, 19352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35852, 0, 3,
                                                                       34817, 18542, 34862, 9092,
                                                                       9128, 19412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35942, 0, 3,
                                                                       34862, 18572, 34907, 9128,
                                                                       9164, 19472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36032, 0, 3,
                                                                       34907, 18602, 34952, 9164,
                                                                       9200, 19532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36122, 0, 3,
                                                                       34997, 18722, 35042, 9272,
                                                                       9308, 19712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36212, 0, 3,
                                                                       35042, 18752, 35087, 9308,
                                                                       9344, 19772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36302, 0, 3,
                                                                       35087, 18782, 35132, 9344,
                                                                       9380, 19832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36392, 0, 3,
                                                                       35132, 18812, 35177, 9380,
                                                                       9416, 19892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36482, 0, 3,
                                                                       35177, 18842, 35222, 9416,
                                                                       9452, 19952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36572, 0, 3,
                                                                       35222, 18872, 35267, 9452,
                                                                       9488, 20012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36662, 0, 3,
                                                                       35267, 18902, 35312, 9488,
                                                                       9524, 20072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36752, 0, 3,
                                                                       35312, 18932, 35357, 9524,
                                                                       9560, 20132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36842, 0, 3,
                                                                       35402, 19112, 35492, 9632,
                                                                       9692, 20392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36992, 0, 3,
                                                                       35492, 19172, 35582, 9692,
                                                                       9752, 20492, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37142, 0, 3,
                                                                       35582, 19232, 35672, 9752,
                                                                       9812, 20592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37292, 0, 3,
                                                                       35672, 19292, 35762, 9812,
                                                                       9872, 20692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37442, 0, 3,
                                                                       35762, 19352, 35852, 9872,
                                                                       9932, 20792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37592, 0, 3,
                                                                       35852, 19412, 35942, 9932,
                                                                       9992, 20892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37742, 0, 3,
                                                                       35942, 19472, 36032, 9992,
                                                                       10052, 20992, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37892, 0, 3,
                                                                       36122, 19712, 36212,
                                                                       10172, 10232, 21292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38042, 0, 3,
                                                                       36212, 19772, 36302,
                                                                       10232, 10292, 21392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38192, 0, 3,
                                                                       36302, 19832, 36392,
                                                                       10292, 10352, 21492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38342, 0, 3,
                                                                       36392, 19892, 36482,
                                                                       10352, 10412, 21592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38492, 0, 3,
                                                                       36482, 19952, 36572,
                                                                       10412, 10472, 21692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38642, 0, 3,
                                                                       36572, 20012, 36662,
                                                                       10472, 10532, 21792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38792, 0, 3,
                                                                       36662, 20072, 36752,
                                                                       10532, 10592, 21892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 38942, 0, 3,
                                                                       36842, 20392, 36992,
                                                                       10712, 10802, 22292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 39167, 0, 3,
                                                                       36992, 20492, 37142,
                                                                       10802, 10892, 22442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 39392, 0, 3,
                                                                       37142, 20592, 37292,
                                                                       10892, 10982, 22592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 39617, 0, 3,
                                                                       37292, 20692, 37442,
                                                                       10982, 11072, 22742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 39842, 0, 3,
                                                                       37442, 20792, 37592,
                                                                       11072, 11162, 22892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40067, 0, 3,
                                                                       37592, 20892, 37742,
                                                                       11162, 11252, 23042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40292, 0, 3,
                                                                       37892, 21292, 38042,
                                                                       11432, 11522, 23492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40517, 0, 3,
                                                                       38042, 21392, 38192,
                                                                       11522, 11612, 23642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40742, 0, 3,
                                                                       38192, 21492, 38342,
                                                                       11612, 11702, 23792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40967, 0, 3,
                                                                       38342, 21592, 38492,
                                                                       11702, 11792, 23942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41192, 0, 3,
                                                                       38492, 21692, 38642,
                                                                       11792, 11882, 24092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41417, 0, 3,
                                                                       38642, 21792, 38792,
                                                                       11882, 11972, 24242,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 41642, 0, 3,
                                                                       38942, 22292, 39167,
                                                                       12152, 12278, 24812,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 41957, 0, 3,
                                                                       39167, 22442, 39392,
                                                                       12278, 12404, 25022,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 42272, 0, 3,
                                                                       39392, 22592, 39617,
                                                                       12404, 12530, 25232,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 42587, 0, 3,
                                                                       39617, 22742, 39842,
                                                                       12530, 12656, 25442,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 42902, 0, 3,
                                                                       39842, 22892, 40067,
                                                                       12656, 12782, 25652,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43217, 0, 3,
                                                                       40292, 23492, 40517,
                                                                       13034, 13160, 26282,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43532, 0, 3,
                                                                       40517, 23642, 40742,
                                                                       13160, 13286, 26492,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43847, 0, 3,
                                                                       40742, 23792, 40967,
                                                                       13286, 13412, 26702,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44162, 0, 3,
                                                                       40967, 23942, 41192,
                                                                       13412, 13538, 26912,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44477, 0, 3,
                                                                       41192, 24092, 41417,
                                                                       13538, 13664, 27122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 44792, 0, 3,
                                                                       41642, 24812, 41957,
                                                                       13916, 14084, 27892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 45212, 0, 3,
                                                                       41957, 25022, 42272,
                                                                       14084, 14252, 28172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 45632, 0, 3,
                                                                       42272, 25232, 42587,
                                                                       14252, 14420, 28452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46052, 0, 3,
                                                                       42587, 25442, 42902,
                                                                       14420, 14588, 28732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46472, 0, 3,
                                                                       43217, 26282, 43532,
                                                                       14924, 15092, 29572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46892, 0, 3,
                                                                       43532, 26492, 43847,
                                                                       15092, 15260, 29852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 47312, 0, 3,
                                                                       43847, 26702, 44162,
                                                                       15260, 15428, 30132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 47732, 0, 3,
                                                                       44162, 26912, 44477,
                                                                       15428, 15596, 30412,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 48152, 0, 3,
                                                                       44792, 27892, 45212,
                                                                       15932, 16148, 31412,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 48692, 0, 3,
                                                                       45212, 28172, 45632,
                                                                       16148, 16364, 31772,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 49232, 0, 3,
                                                                       45632, 28452, 46052,
                                                                       16364, 16580, 32132,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 49772, 0, 3,
                                                                       46472, 29572, 46892,
                                                                       17012, 17228, 33212,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 50312, 0, 3,
                                                                       46892, 29852, 47312,
                                                                       17228, 17444, 33572,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 50852, 0, 3,
                                                                       47312, 30132, 47732,
                                                                       17444, 17660, 33932,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51392, 3, 18092,
                                                                       18102, 34292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51413, 3, 18102,
                                                                       18112, 34307, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51434, 3, 18112,
                                                                       18122, 34322, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51455, 3, 18122,
                                                                       18132, 34337, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51476, 3, 18132,
                                                                       18142, 34352, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51497, 3, 18142,
                                                                       18152, 34367, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51518, 3, 18152,
                                                                       18162, 34382, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51539, 3, 18162,
                                                                       18172, 34397, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51560, 3, 18172,
                                                                       18182, 34412, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51581, 3, 18182,
                                                                       18192, 34427, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51602, 3, 18212,
                                                                       18222, 34442, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51623, 3, 18222,
                                                                       18232, 34457, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51644, 3, 18232,
                                                                       18242, 34472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51665, 3, 18242,
                                                                       18252, 34487, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51686, 3, 18252,
                                                                       18262, 34502, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51707, 3, 18262,
                                                                       18272, 34517, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51728, 3, 18272,
                                                                       18282, 34532, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51749, 3, 18282,
                                                                       18292, 34547, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51770, 3, 18292,
                                                                       18302, 34562, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51791, 3, 18302,
                                                                       18312, 34577, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51812, 0, 3,
                                                                       51392, 34292, 51413,
                                                                       18332, 18362, 34592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51875, 0, 3,
                                                                       51413, 34307, 51434,
                                                                       18362, 18392, 34637,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51938, 0, 3,
                                                                       51434, 34322, 51455,
                                                                       18392, 18422, 34682,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52001, 0, 3,
                                                                       51455, 34337, 51476,
                                                                       18422, 18452, 34727,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52064, 0, 3,
                                                                       51476, 34352, 51497,
                                                                       18452, 18482, 34772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52127, 0, 3,
                                                                       51497, 34367, 51518,
                                                                       18482, 18512, 34817,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52190, 0, 3,
                                                                       51518, 34382, 51539,
                                                                       18512, 18542, 34862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52253, 0, 3,
                                                                       51539, 34397, 51560,
                                                                       18542, 18572, 34907,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52316, 0, 3,
                                                                       51560, 34412, 51581,
                                                                       18572, 18602, 34952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52379, 0, 3,
                                                                       51602, 34442, 51623,
                                                                       18662, 18692, 34997,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52442, 0, 3,
                                                                       51623, 34457, 51644,
                                                                       18692, 18722, 35042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52505, 0, 3,
                                                                       51644, 34472, 51665,
                                                                       18722, 18752, 35087,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52568, 0, 3,
                                                                       51665, 34487, 51686,
                                                                       18752, 18782, 35132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52631, 0, 3,
                                                                       51686, 34502, 51707,
                                                                       18782, 18812, 35177,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52694, 0, 3,
                                                                       51707, 34517, 51728,
                                                                       18812, 18842, 35222,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52757, 0, 3,
                                                                       51728, 34532, 51749,
                                                                       18842, 18872, 35267,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52820, 0, 3,
                                                                       51749, 34547, 51770,
                                                                       18872, 18902, 35312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 52883, 0, 3,
                                                                       51770, 34562, 51791,
                                                                       18902, 18932, 35357,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52946, 0, 3,
                                                                       51812, 34592, 51875,
                                                                       18992, 19052, 35402,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53072, 0, 3,
                                                                       51875, 34637, 51938,
                                                                       19052, 19112, 35492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53198, 0, 3,
                                                                       51938, 34682, 52001,
                                                                       19112, 19172, 35582,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53324, 0, 3,
                                                                       52001, 34727, 52064,
                                                                       19172, 19232, 35672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53450, 0, 3,
                                                                       52064, 34772, 52127,
                                                                       19232, 19292, 35762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53576, 0, 3,
                                                                       52127, 34817, 52190,
                                                                       19292, 19352, 35852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53702, 0, 3,
                                                                       52190, 34862, 52253,
                                                                       19352, 19412, 35942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53828, 0, 3,
                                                                       52253, 34907, 52316,
                                                                       19412, 19472, 36032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53954, 0, 3,
                                                                       52379, 34997, 52442,
                                                                       19592, 19652, 36122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 54080, 0, 3,
                                                                       52442, 35042, 52505,
                                                                       19652, 19712, 36212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 54206, 0, 3,
                                                                       52505, 35087, 52568,
                                                                       19712, 19772, 36302,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 54332, 0, 3,
                                                                       52568, 35132, 52631,
                                                                       19772, 19832, 36392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 54458, 0, 3,
                                                                       52631, 35177, 52694,
                                                                       19832, 19892, 36482,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 54584, 0, 3,
                                                                       52694, 35222, 52757,
                                                                       19892, 19952, 36572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 54710, 0, 3,
                                                                       52757, 35267, 52820,
                                                                       19952, 20012, 36662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 54836, 0, 3,
                                                                       52820, 35312, 52883,
                                                                       20012, 20072, 36752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54962, 0, 3,
                                                                       52946, 35402, 53072,
                                                                       20192, 20292, 36842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 55172, 0, 3,
                                                                       53072, 35492, 53198,
                                                                       20292, 20392, 36992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 55382, 0, 3,
                                                                       53198, 35582, 53324,
                                                                       20392, 20492, 37142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 55592, 0, 3,
                                                                       53324, 35672, 53450,
                                                                       20492, 20592, 37292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 55802, 0, 3,
                                                                       53450, 35762, 53576,
                                                                       20592, 20692, 37442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 56012, 0, 3,
                                                                       53576, 35852, 53702,
                                                                       20692, 20792, 37592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 56222, 0, 3,
                                                                       53702, 35942, 53828,
                                                                       20792, 20892, 37742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 56432, 0, 3,
                                                                       53954, 36122, 54080,
                                                                       21092, 21192, 37892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 56642, 0, 3,
                                                                       54080, 36212, 54206,
                                                                       21192, 21292, 38042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 56852, 0, 3,
                                                                       54206, 36302, 54332,
                                                                       21292, 21392, 38192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 57062, 0, 3,
                                                                       54332, 36392, 54458,
                                                                       21392, 21492, 38342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 57272, 0, 3,
                                                                       54458, 36482, 54584,
                                                                       21492, 21592, 38492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 57482, 0, 3,
                                                                       54584, 36572, 54710,
                                                                       21592, 21692, 38642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 57692, 0, 3,
                                                                       54710, 36662, 54836,
                                                                       21692, 21792, 38792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 57902, 0, 3,
                                                                       54962, 36842, 55172,
                                                                       21992, 22142, 38942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 58217, 0, 3,
                                                                       55172, 36992, 55382,
                                                                       22142, 22292, 39167,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 58532, 0, 3,
                                                                       55382, 37142, 55592,
                                                                       22292, 22442, 39392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 58847, 0, 3,
                                                                       55592, 37292, 55802,
                                                                       22442, 22592, 39617,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 59162, 0, 3,
                                                                       55802, 37442, 56012,
                                                                       22592, 22742, 39842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 59477, 0, 3,
                                                                       56012, 37592, 56222,
                                                                       22742, 22892, 40067,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 59792, 0, 3,
                                                                       56432, 37892, 56642,
                                                                       23192, 23342, 40292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 60107, 0, 3,
                                                                       56642, 38042, 56852,
                                                                       23342, 23492, 40517,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 60422, 0, 3,
                                                                       56852, 38192, 57062,
                                                                       23492, 23642, 40742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 60737, 0, 3,
                                                                       57062, 38342, 57272,
                                                                       23642, 23792, 40967,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 61052, 0, 3,
                                                                       57272, 38492, 57482,
                                                                       23792, 23942, 41192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 61367, 0, 3,
                                                                       57482, 38642, 57692,
                                                                       23942, 24092, 41417,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 61682, 0, 3,
                                                                       57902, 38942, 58217,
                                                                       24392, 24602, 41642,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 62123, 0, 3,
                                                                       58217, 39167, 58532,
                                                                       24602, 24812, 41957,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 62564, 0, 3,
                                                                       58532, 39392, 58847,
                                                                       24812, 25022, 42272,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 63005, 0, 3,
                                                                       58847, 39617, 59162,
                                                                       25022, 25232, 42587,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 63446, 0, 3,
                                                                       59162, 39842, 59477,
                                                                       25232, 25442, 42902,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 63887, 0, 3,
                                                                       59792, 40292, 60107,
                                                                       25862, 26072, 43217,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 64328, 0, 3,
                                                                       60107, 40517, 60422,
                                                                       26072, 26282, 43532,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 64769, 0, 3,
                                                                       60422, 40742, 60737,
                                                                       26282, 26492, 43847,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 65210, 0, 3,
                                                                       60737, 40967, 61052,
                                                                       26492, 26702, 44162,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 65651, 0, 3,
                                                                       61052, 41192, 61367,
                                                                       26702, 26912, 44477,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 66092, 0, 3,
                                                                       61682, 41642, 62123,
                                                                       27332, 27612, 44792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 66680, 0, 3,
                                                                       62123, 41957, 62564,
                                                                       27612, 27892, 45212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 67268, 0, 3,
                                                                       62564, 42272, 63005,
                                                                       27892, 28172, 45632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 67856, 0, 3,
                                                                       63005, 42587, 63446,
                                                                       28172, 28452, 46052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 68444, 0, 3,
                                                                       63887, 43217, 64328,
                                                                       29012, 29292, 46472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 69032, 0, 3,
                                                                       64328, 43532, 64769,
                                                                       29292, 29572, 46892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 69620, 0, 3,
                                                                       64769, 43847, 65210,
                                                                       29572, 29852, 47312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 70208, 0, 3,
                                                                       65210, 44162, 65651,
                                                                       29852, 30132, 47732,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 70796, 0, 3,
                                                                       66092, 44792, 66680,
                                                                       30692, 31052, 48152,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 71552, 0, 3,
                                                                       66680, 45212, 67268,
                                                                       31052, 31412, 48692,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 72308, 0, 3,
                                                                       67268, 45632, 67856,
                                                                       31412, 31772, 49232,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 73064, 0, 3,
                                                                       68444, 46472, 69032,
                                                                       32492, 32852, 49772,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 73820, 0, 3,
                                                                       69032, 46892, 69620,
                                                                       32852, 33212, 50312,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 74576, 0, 3,
                                                                       69620, 47312, 70208,
                                                                       33212, 33572, 50852,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75332, 3, 34292,
                                                                       34307, 51434, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75360, 3, 34307,
                                                                       34322, 51455, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75388, 3, 34322,
                                                                       34337, 51476, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75416, 3, 34337,
                                                                       34352, 51497, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75444, 3, 34352,
                                                                       34367, 51518, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75472, 3, 34367,
                                                                       34382, 51539, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75500, 3, 34382,
                                                                       34397, 51560, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75528, 3, 34397,
                                                                       34412, 51581, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75556, 3, 34442,
                                                                       34457, 51644, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75584, 3, 34457,
                                                                       34472, 51665, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75612, 3, 34472,
                                                                       34487, 51686, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75640, 3, 34487,
                                                                       34502, 51707, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75668, 3, 34502,
                                                                       34517, 51728, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75696, 3, 34517,
                                                                       34532, 51749, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75724, 3, 34532,
                                                                       34547, 51770, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75752, 3, 34547,
                                                                       34562, 51791, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75780, 0, 3,
                                                                       75332, 51434, 75360,
                                                                       34592, 34637, 51938,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75864, 0, 3,
                                                                       75360, 51455, 75388,
                                                                       34637, 34682, 52001,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75948, 0, 3,
                                                                       75388, 51476, 75416,
                                                                       34682, 34727, 52064,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76032, 0, 3,
                                                                       75416, 51497, 75444,
                                                                       34727, 34772, 52127,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76116, 0, 3,
                                                                       75444, 51518, 75472,
                                                                       34772, 34817, 52190,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76200, 0, 3,
                                                                       75472, 51539, 75500,
                                                                       34817, 34862, 52253,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76284, 0, 3,
                                                                       75500, 51560, 75528,
                                                                       34862, 34907, 52316,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76368, 0, 3,
                                                                       75556, 51644, 75584,
                                                                       34997, 35042, 52505,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76452, 0, 3,
                                                                       75584, 51665, 75612,
                                                                       35042, 35087, 52568,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76536, 0, 3,
                                                                       75612, 51686, 75640,
                                                                       35087, 35132, 52631,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76620, 0, 3,
                                                                       75640, 51707, 75668,
                                                                       35132, 35177, 52694,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76704, 0, 3,
                                                                       75668, 51728, 75696,
                                                                       35177, 35222, 52757,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76788, 0, 3,
                                                                       75696, 51749, 75724,
                                                                       35222, 35267, 52820,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76872, 0, 3,
                                                                       75724, 51770, 75752,
                                                                       35267, 35312, 52883,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76956, 0, 3,
                                                                       75780, 51938, 75864,
                                                                       35402, 35492, 53198,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77124, 0, 3,
                                                                       75864, 52001, 75948,
                                                                       35492, 35582, 53324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77292, 0, 3,
                                                                       75948, 52064, 76032,
                                                                       35582, 35672, 53450,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77460, 0, 3,
                                                                       76032, 52127, 76116,
                                                                       35672, 35762, 53576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77628, 0, 3,
                                                                       76116, 52190, 76200,
                                                                       35762, 35852, 53702,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77796, 0, 3,
                                                                       76200, 52253, 76284,
                                                                       35852, 35942, 53828,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77964, 0, 3,
                                                                       76368, 52505, 76452,
                                                                       36122, 36212, 54206,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 78132, 0, 3,
                                                                       76452, 52568, 76536,
                                                                       36212, 36302, 54332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 78300, 0, 3,
                                                                       76536, 52631, 76620,
                                                                       36302, 36392, 54458,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 78468, 0, 3,
                                                                       76620, 52694, 76704,
                                                                       36392, 36482, 54584,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 78636, 0, 3,
                                                                       76704, 52757, 76788,
                                                                       36482, 36572, 54710,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 78804, 0, 3,
                                                                       76788, 52820, 76872,
                                                                       36572, 36662, 54836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 78972, 0, 3,
                                                                       76956, 53198, 77124,
                                                                       36842, 36992, 55382,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 79252, 0, 3,
                                                                       77124, 53324, 77292,
                                                                       36992, 37142, 55592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 79532, 0, 3,
                                                                       77292, 53450, 77460,
                                                                       37142, 37292, 55802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 79812, 0, 3,
                                                                       77460, 53576, 77628,
                                                                       37292, 37442, 56012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 80092, 0, 3,
                                                                       77628, 53702, 77796,
                                                                       37442, 37592, 56222,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 80372, 0, 3,
                                                                       77964, 54206, 78132,
                                                                       37892, 38042, 56852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 80652, 0, 3,
                                                                       78132, 54332, 78300,
                                                                       38042, 38192, 57062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 80932, 0, 3,
                                                                       78300, 54458, 78468,
                                                                       38192, 38342, 57272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 81212, 0, 3,
                                                                       78468, 54584, 78636,
                                                                       38342, 38492, 57482,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 81492, 0, 3,
                                                                       78636, 54710, 78804,
                                                                       38492, 38642, 57692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 81772, 0, 3,
                                                                       78972, 55382, 79252,
                                                                       38942, 39167, 58532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 82192, 0, 3,
                                                                       79252, 55592, 79532,
                                                                       39167, 39392, 58847,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 82612, 0, 3,
                                                                       79532, 55802, 79812,
                                                                       39392, 39617, 59162,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 83032, 0, 3,
                                                                       79812, 56012, 80092,
                                                                       39617, 39842, 59477,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 83452, 0, 3,
                                                                       80372, 56852, 80652,
                                                                       40292, 40517, 60422,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 83872, 0, 3,
                                                                       80652, 57062, 80932,
                                                                       40517, 40742, 60737,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 84292, 0, 3,
                                                                       80932, 57272, 81212,
                                                                       40742, 40967, 61052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 84712, 0, 3,
                                                                       81212, 57482, 81492,
                                                                       40967, 41192, 61367,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 85132, 0, 3,
                                                                       81772, 58532, 82192,
                                                                       41642, 41957, 62564,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 85720, 0, 3,
                                                                       82192, 58847, 82612,
                                                                       41957, 42272, 63005,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 86308, 0, 3,
                                                                       82612, 59162, 83032,
                                                                       42272, 42587, 63446,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 86896, 0, 3,
                                                                       83452, 60422, 83872,
                                                                       43217, 43532, 64769,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 87484, 0, 3,
                                                                       83872, 60737, 84292,
                                                                       43532, 43847, 65210,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 88072, 0, 3,
                                                                       84292, 61052, 84712,
                                                                       43847, 44162, 65651,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 88660, 0, 3,
                                                                       85132, 62564, 85720,
                                                                       44792, 45212, 67268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 89444, 0, 3,
                                                                       85720, 63005, 86308,
                                                                       45212, 45632, 67856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 90228, 0, 3,
                                                                       86896, 64769, 87484,
                                                                       46472, 46892, 69620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 91012, 0, 3,
                                                                       87484, 65210, 88072,
                                                                       46892, 47312, 70208,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 91796, 0, 3,
                                                                       88660, 67268, 89444,
                                                                       48152, 48692, 72308,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 92804, 0, 3,
                                                                       90228, 69620, 91012,
                                                                       49772, 50312, 74576,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93812, 3, 51392,
                                                                       51413, 75332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93848, 3, 51413,
                                                                       51434, 75360, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93884, 3, 51434,
                                                                       51455, 75388, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93920, 3, 51455,
                                                                       51476, 75416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93956, 3, 51476,
                                                                       51497, 75444, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 93992, 3, 51497,
                                                                       51518, 75472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94028, 3, 51518,
                                                                       51539, 75500, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94064, 3, 51539,
                                                                       51560, 75528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94100, 3, 51602,
                                                                       51623, 75556, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94136, 3, 51623,
                                                                       51644, 75584, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94172, 3, 51644,
                                                                       51665, 75612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94208, 3, 51665,
                                                                       51686, 75640, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94244, 3, 51686,
                                                                       51707, 75668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94280, 3, 51707,
                                                                       51728, 75696, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94316, 3, 51728,
                                                                       51749, 75724, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 94352, 3, 51749,
                                                                       51770, 75752, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 94388, 0, 3,
                                                                       93812, 75332, 93848,
                                                                       51812, 51875, 75780,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 94496, 0, 3,
                                                                       93848, 75360, 93884,
                                                                       51875, 51938, 75864,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 94604, 0, 3,
                                                                       93884, 75388, 93920,
                                                                       51938, 52001, 75948,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 94712, 0, 3,
                                                                       93920, 75416, 93956,
                                                                       52001, 52064, 76032,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 94820, 0, 3,
                                                                       93956, 75444, 93992,
                                                                       52064, 52127, 76116,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 94928, 0, 3,
                                                                       93992, 75472, 94028,
                                                                       52127, 52190, 76200,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95036, 0, 3,
                                                                       94028, 75500, 94064,
                                                                       52190, 52253, 76284,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95144, 0, 3,
                                                                       94100, 75556, 94136,
                                                                       52379, 52442, 76368,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95252, 0, 3,
                                                                       94136, 75584, 94172,
                                                                       52442, 52505, 76452,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95360, 0, 3,
                                                                       94172, 75612, 94208,
                                                                       52505, 52568, 76536,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95468, 0, 3,
                                                                       94208, 75640, 94244,
                                                                       52568, 52631, 76620,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95576, 0, 3,
                                                                       94244, 75668, 94280,
                                                                       52631, 52694, 76704,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95684, 0, 3,
                                                                       94280, 75696, 94316,
                                                                       52694, 52757, 76788,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95792, 0, 3,
                                                                       94316, 75724, 94352,
                                                                       52757, 52820, 76872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 95900, 0, 3,
                                                                       94388, 75780, 94496,
                                                                       52946, 53072, 76956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 96116, 0, 3,
                                                                       94496, 75864, 94604,
                                                                       53072, 53198, 77124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 96332, 0, 3,
                                                                       94604, 75948, 94712,
                                                                       53198, 53324, 77292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 96548, 0, 3,
                                                                       94712, 76032, 94820,
                                                                       53324, 53450, 77460,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 96764, 0, 3,
                                                                       94820, 76116, 94928,
                                                                       53450, 53576, 77628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 96980, 0, 3,
                                                                       94928, 76200, 95036,
                                                                       53576, 53702, 77796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97196, 0, 3,
                                                                       95144, 76368, 95252,
                                                                       53954, 54080, 77964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97412, 0, 3,
                                                                       95252, 76452, 95360,
                                                                       54080, 54206, 78132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97628, 0, 3,
                                                                       95360, 76536, 95468,
                                                                       54206, 54332, 78300,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97844, 0, 3,
                                                                       95468, 76620, 95576,
                                                                       54332, 54458, 78468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 98060, 0, 3,
                                                                       95576, 76704, 95684,
                                                                       54458, 54584, 78636,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 98276, 0, 3,
                                                                       95684, 76788, 95792,
                                                                       54584, 54710, 78804,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 98492, 0, 3,
                                                                       95900, 76956, 96116,
                                                                       54962, 55172, 78972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 98852, 0, 3,
                                                                       96116, 77124, 96332,
                                                                       55172, 55382, 79252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 99212, 0, 3,
                                                                       96332, 77292, 96548,
                                                                       55382, 55592, 79532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 99572, 0, 3,
                                                                       96548, 77460, 96764,
                                                                       55592, 55802, 79812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 99932, 0, 3,
                                                                       96764, 77628, 96980,
                                                                       55802, 56012, 80092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 100292, 0, 3,
                                                                       97196, 77964, 97412,
                                                                       56432, 56642, 80372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 100652, 0, 3,
                                                                       97412, 78132, 97628,
                                                                       56642, 56852, 80652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 101012, 0, 3,
                                                                       97628, 78300, 97844,
                                                                       56852, 57062, 80932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 101372, 0, 3,
                                                                       97844, 78468, 98060,
                                                                       57062, 57272, 81212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 101732, 0, 3,
                                                                       98060, 78636, 98276,
                                                                       57272, 57482, 81492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 102092, 0, 3,
                                                                       98492, 78972, 98852,
                                                                       57902, 58217, 81772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 102632, 0, 3,
                                                                       98852, 79252, 99212,
                                                                       58217, 58532, 82192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 103172, 0, 3,
                                                                       99212, 79532, 99572,
                                                                       58532, 58847, 82612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 103712, 0, 3,
                                                                       99572, 79812, 99932,
                                                                       58847, 59162, 83032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 104252, 0, 3,
                                                                       100292, 80372, 100652,
                                                                       59792, 60107, 83452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 104792, 0, 3,
                                                                       100652, 80652, 101012,
                                                                       60107, 60422, 83872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 105332, 0, 3,
                                                                       101012, 80932, 101372,
                                                                       60422, 60737, 84292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 105872, 0, 3,
                                                                       101372, 81212, 101732,
                                                                       60737, 61052, 84712,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 106412, 0, 3,
                                                                       102092, 81772, 102632,
                                                                       61682, 62123, 85132,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 107168, 0, 3,
                                                                       102632, 82192, 103172,
                                                                       62123, 62564, 85720,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 107924, 0, 3,
                                                                       103172, 82612, 103712,
                                                                       62564, 63005, 86308,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 108680, 0, 3,
                                                                       104252, 83452, 104792,
                                                                       63887, 64328, 86896,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 109436, 0, 3,
                                                                       104792, 83872, 105332,
                                                                       64328, 64769, 87484,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 110192, 0, 3,
                                                                       105332, 84292, 105872,
                                                                       64769, 65210, 88072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 110948, 0, 3,
                                                                       106412, 85132, 107168,
                                                                       66092, 66680, 88660,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 111956, 0, 3,
                                                                       107168, 85720, 107924,
                                                                       66680, 67268, 89444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 112964, 0, 3,
                                                                       108680, 86896, 109436,
                                                                       68444, 69032, 90228,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 113972, 0, 3,
                                                                       109436, 87484, 110192,
                                                                       69032, 69620, 91012,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 114980, 0, 3,
                                                                       110948, 88660, 111956,
                                                                       70796, 71552, 91796,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 116276, 0, 3,
                                                                       112964, 90228, 113972,
                                                                       73064, 73820, 92804,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 117572, 102092, 540, ncols);

                    simdfunc::contract_primitives(buffer, 118337, 104252, 540, ncols);

                    simdfunc::contract_primitives(buffer, 119102, 106412, 756, ncols);

                    simdfunc::contract_primitives(buffer, 120173, 108680, 756, ncols);

                    simdfunc::contract_primitives(buffer, 121244, 110948, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 122672, 112964, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 124100, 114980, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 125936, 116276, 1296, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 118112, 117572, 15, 1, nmax);

        simdtrf::transform_k_inner(buffer, 118877, 118337, 15, 1, nmax);

        simdtrf::transform_k_inner(buffer, 119858, 119102, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 120929, 120173, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 122252, 121244, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 123680, 122672, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 125396, 124100, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 127232, 125936, 36, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 127772, 118112, 119858, 15, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 128447, 118877, 120929, 15, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 129122, 119858, 122252, 15, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 130067, 120929, 123680, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 131012, 122252, 125396, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 132272, 123680, 127232, 15, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 133532, 127772, 129122, 15, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 134882, 128447, 130067, 15, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 136232, 129122, 131012, 15, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 138122, 130067, 132272, 15, nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 140012, 133532, 136232, 15, nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 142262, 134882, 138122, 15, nmax);

        simdtrf::transform_g_inner(buffer, 144512, 142262, 10, 15, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 144512, 135, nmax);

        simdtrf::transform_g_inner(buffer, 144512, 140012, 10, 15, nmax);

        simdtrf::transform_f_outer(values + 945 * nvalues + n * npairs, nvalues, buffer, 144512,
                                   135, nmax);
    }

    for (size_t m = 0; m < 1890; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
