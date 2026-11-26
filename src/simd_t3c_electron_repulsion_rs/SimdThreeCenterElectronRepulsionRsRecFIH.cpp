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


#include "SimdThreeCenterElectronRepulsionRsRecFIH.hpp"

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
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_fih_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_fih_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 134900, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 134900, 101172, 9891, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3198, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3201, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3204, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3207, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3210, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3213, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3216, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3219, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3222, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3225, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3228, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3231, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3234, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3237, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3240, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3243, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3246, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3249, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3252, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3255, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3258, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3261, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3264, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3267, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3270, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3273, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3276, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3279, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3282, 3, 9, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3291, 3, 10, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3300, 3, 11, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3309, 3, 12, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3318, 3, 13, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3327, 3, 14, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3336, 3, 15, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3345, 3, 16, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3354, 3, 17, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3363, 3, 18, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3372, 3, 19, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3381, 3, 24, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3390, 3, 25, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3399, 3, 26, 87,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3408, 3, 27, 90,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3417, 3, 28, 93,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3426, 3, 29, 96,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3435, 3, 30, 99,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3444, 3, 31, 102,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3453, 3, 32, 105,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3462, 3, 33, 108,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3471, 3, 34, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3480, 3, 36, 114,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3498, 3, 39, 120,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3516, 3, 42, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3534, 3, 45, 132,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3552, 3, 48, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3570, 3, 51, 144,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3588, 3, 54, 150,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3606, 3, 57, 156,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3624, 3, 60, 162,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3642, 3, 63, 168,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3660, 3, 66, 174,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3678, 3, 69, 180,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3696, 3, 75, 186,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3714, 3, 78, 192,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3732, 3, 81, 198,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3750, 3, 84, 204,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3768, 3, 87, 210,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3786, 3, 90, 216,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3804, 3, 93, 222,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3822, 3, 96, 228,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3840, 3, 99, 234,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3858, 3, 102, 240,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3876, 3, 105, 246,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3894, 3, 108, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3912, 3, 114, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3942, 3, 120, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3972, 3, 126, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4002, 3, 132, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4032, 3, 138, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4062, 3, 144, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4092, 3, 150, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4122, 3, 156, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4152, 3, 162, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4182, 3, 168, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4212, 3, 174, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4242, 3, 186, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4272, 3, 192, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4302, 3, 198, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4332, 3, 204, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4362, 3, 210, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4392, 3, 216, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4422, 3, 222, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4452, 3, 228, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4482, 3, 234, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4512, 3, 240, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4542, 3, 246, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4572, 3, 258, 478,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4617, 3, 268, 493,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4662, 3, 278, 508,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4707, 3, 288, 523,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4752, 3, 298, 538,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4797, 3, 308, 553,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4842, 3, 318, 568,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4887, 3, 328, 583,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4932, 3, 338, 598,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4977, 3, 348, 613,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5022, 3, 368, 628,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5067, 3, 378, 643,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5112, 3, 388, 658,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5157, 3, 398, 673,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5202, 3, 408, 688,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5247, 3, 418, 703,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5292, 3, 428, 718,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5337, 3, 438, 733,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5382, 3, 448, 748,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5427, 3, 458, 763,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5472, 3, 478, 778,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5535, 3, 493, 799,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5598, 3, 508, 820,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5661, 3, 523, 841,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5724, 3, 538, 862,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5787, 3, 553, 883,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5850, 3, 568, 904,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5913, 3, 583, 925,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5976, 3, 598, 946,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6039, 3, 628, 967,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6102, 3, 643, 988,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6165, 3, 658,
                                                                       1009, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6228, 3, 673,
                                                                       1030, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6291, 3, 688,
                                                                       1051, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6354, 3, 703,
                                                                       1072, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6417, 3, 718,
                                                                       1093, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6480, 3, 733,
                                                                       1114, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6543, 3, 748,
                                                                       1135, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6606, 3, 778,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6690, 3, 799,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6774, 3, 820,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6858, 3, 841,
                                                                       1240, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6942, 3, 862,
                                                                       1268, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7026, 3, 883,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7110, 3, 904,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7194, 3, 925,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7278, 3, 967,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7362, 3, 988,
                                                                       1408, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7446, 3, 1009,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7530, 3, 1030,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7614, 3, 1051,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7698, 3, 1072,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7782, 3, 1093,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7866, 3, 1114,
                                                                       1576, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7950, 3, 1156,
                                                                       1604, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8058, 3, 1184,
                                                                       1640, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8166, 3, 1212,
                                                                       1676, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8274, 3, 1240,
                                                                       1712, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8382, 3, 1268,
                                                                       1748, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8490, 3, 1296,
                                                                       1784, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8598, 3, 1324,
                                                                       1820, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8706, 3, 1380,
                                                                       1856, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8814, 3, 1408,
                                                                       1892, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8922, 3, 1436,
                                                                       1928, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9030, 3, 1464,
                                                                       1964, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9138, 3, 1492,
                                                                       2000, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9246, 3, 1520,
                                                                       2036, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9354, 3, 1548,
                                                                       2072, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9462, 3, 1604,
                                                                       2108, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9597, 3, 1640,
                                                                       2153, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9732, 3, 1676,
                                                                       2198, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9867, 3, 1712,
                                                                       2243, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10002, 3, 1748,
                                                                       2288, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10137, 3, 1784,
                                                                       2333, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10272, 3, 1856,
                                                                       2378, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10407, 3, 1892,
                                                                       2423, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10542, 3, 1928,
                                                                       2468, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10677, 3, 1964,
                                                                       2513, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10812, 3, 2000,
                                                                       2558, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10947, 3, 2036,
                                                                       2603, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11082, 3, 2108,
                                                                       2648, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11247, 3, 2153,
                                                                       2703, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11412, 3, 2198,
                                                                       2758, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11577, 3, 2243,
                                                                       2813, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11742, 3, 2288,
                                                                       2868, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11907, 3, 2378,
                                                                       2923, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12072, 3, 2423,
                                                                       2978, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12237, 3, 2468,
                                                                       3033, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12402, 3, 2513,
                                                                       3088, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 12567, 3, 2558,
                                                                       3143, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12732, 3, 7, 8,
                                                                       3204, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12738, 3, 8, 9,
                                                                       3207, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12744, 3, 9, 10,
                                                                       3210, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12750, 3, 10, 11,
                                                                       3213, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12756, 3, 11, 12,
                                                                       3216, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12762, 3, 12, 13,
                                                                       3219, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12768, 3, 13, 14,
                                                                       3222, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12774, 3, 14, 15,
                                                                       3225, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12780, 3, 15, 16,
                                                                       3228, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12786, 3, 16, 17,
                                                                       3231, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12792, 3, 17, 18,
                                                                       3234, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12798, 3, 18, 19,
                                                                       3237, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12804, 3, 22, 23,
                                                                       3246, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12810, 3, 23, 24,
                                                                       3249, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12816, 3, 24, 25,
                                                                       3252, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12822, 3, 25, 26,
                                                                       3255, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12828, 3, 26, 27,
                                                                       3258, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12834, 3, 27, 28,
                                                                       3261, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12840, 3, 28, 29,
                                                                       3264, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12846, 3, 29, 30,
                                                                       3267, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12852, 3, 30, 31,
                                                                       3270, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12858, 3, 31, 32,
                                                                       3273, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12864, 3, 32, 33,
                                                                       3276, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12870, 3, 33, 34,
                                                                       3279, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12876, 0, 3,
                                                                       12732, 3204, 12738, 3282,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12894, 0, 3,
                                                                       12738, 3207, 12744, 3291,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12912, 0, 3,
                                                                       12744, 3210, 12750, 3300,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12930, 0, 3,
                                                                       12750, 3213, 12756, 3309,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12948, 0, 3,
                                                                       12756, 3216, 12762, 3318,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12966, 0, 3,
                                                                       12762, 3219, 12768, 3327,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12984, 0, 3,
                                                                       12768, 3222, 12774, 3336,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13002, 0, 3,
                                                                       12774, 3225, 12780, 3345,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13020, 0, 3,
                                                                       12780, 3228, 12786, 3354,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13038, 0, 3,
                                                                       12786, 3231, 12792, 3363,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13056, 0, 3,
                                                                       12792, 3234, 12798, 3372,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13074, 0, 3,
                                                                       12804, 3246, 12810, 3381,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13092, 0, 3,
                                                                       12810, 3249, 12816, 3390,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13110, 0, 3,
                                                                       12816, 3252, 12822, 3399,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13128, 0, 3,
                                                                       12822, 3255, 12828, 3408,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13146, 0, 3,
                                                                       12828, 3258, 12834, 3417,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13164, 0, 3,
                                                                       12834, 3261, 12840, 3426,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13182, 0, 3,
                                                                       12840, 3264, 12846, 3435,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13200, 0, 3,
                                                                       12846, 3267, 12852, 3444,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13218, 0, 3,
                                                                       12852, 3270, 12858, 3453,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13236, 0, 3,
                                                                       12858, 3273, 12864, 3462,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13254, 0, 3,
                                                                       12864, 3276, 12870, 3471,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13272, 0, 3,
                                                                       12876, 3282, 12894, 114,
                                                                       120, 3516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13308, 0, 3,
                                                                       12894, 3291, 12912, 120,
                                                                       126, 3534, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13344, 0, 3,
                                                                       12912, 3300, 12930, 126,
                                                                       132, 3552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13380, 0, 3,
                                                                       12930, 3309, 12948, 132,
                                                                       138, 3570, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13416, 0, 3,
                                                                       12948, 3318, 12966, 138,
                                                                       144, 3588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13452, 0, 3,
                                                                       12966, 3327, 12984, 144,
                                                                       150, 3606, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13488, 0, 3,
                                                                       12984, 3336, 13002, 150,
                                                                       156, 3624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13524, 0, 3,
                                                                       13002, 3345, 13020, 156,
                                                                       162, 3642, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13560, 0, 3,
                                                                       13020, 3354, 13038, 162,
                                                                       168, 3660, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13596, 0, 3,
                                                                       13038, 3363, 13056, 168,
                                                                       174, 3678, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13632, 0, 3,
                                                                       13074, 3381, 13092, 186,
                                                                       192, 3732, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13668, 0, 3,
                                                                       13092, 3390, 13110, 192,
                                                                       198, 3750, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13704, 0, 3,
                                                                       13110, 3399, 13128, 198,
                                                                       204, 3768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13740, 0, 3,
                                                                       13128, 3408, 13146, 204,
                                                                       210, 3786, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13776, 0, 3,
                                                                       13146, 3417, 13164, 210,
                                                                       216, 3804, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13812, 0, 3,
                                                                       13164, 3426, 13182, 216,
                                                                       222, 3822, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13848, 0, 3,
                                                                       13182, 3435, 13200, 222,
                                                                       228, 3840, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13884, 0, 3,
                                                                       13200, 3444, 13218, 228,
                                                                       234, 3858, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13920, 0, 3,
                                                                       13218, 3453, 13236, 234,
                                                                       240, 3876, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13956, 0, 3,
                                                                       13236, 3462, 13254, 240,
                                                                       246, 3894, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13992, 0, 3,
                                                                       13272, 3516, 13308, 258,
                                                                       268, 3972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14052, 0, 3,
                                                                       13308, 3534, 13344, 268,
                                                                       278, 4002, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14112, 0, 3,
                                                                       13344, 3552, 13380, 278,
                                                                       288, 4032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14172, 0, 3,
                                                                       13380, 3570, 13416, 288,
                                                                       298, 4062, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14232, 0, 3,
                                                                       13416, 3588, 13452, 298,
                                                                       308, 4092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14292, 0, 3,
                                                                       13452, 3606, 13488, 308,
                                                                       318, 4122, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14352, 0, 3,
                                                                       13488, 3624, 13524, 318,
                                                                       328, 4152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14412, 0, 3,
                                                                       13524, 3642, 13560, 328,
                                                                       338, 4182, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14472, 0, 3,
                                                                       13560, 3660, 13596, 338,
                                                                       348, 4212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14532, 0, 3,
                                                                       13632, 3732, 13668, 368,
                                                                       378, 4302, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14592, 0, 3,
                                                                       13668, 3750, 13704, 378,
                                                                       388, 4332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14652, 0, 3,
                                                                       13704, 3768, 13740, 388,
                                                                       398, 4362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14712, 0, 3,
                                                                       13740, 3786, 13776, 398,
                                                                       408, 4392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14772, 0, 3,
                                                                       13776, 3804, 13812, 408,
                                                                       418, 4422, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14832, 0, 3,
                                                                       13812, 3822, 13848, 418,
                                                                       428, 4452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14892, 0, 3,
                                                                       13848, 3840, 13884, 428,
                                                                       438, 4482, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14952, 0, 3,
                                                                       13884, 3858, 13920, 438,
                                                                       448, 4512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15012, 0, 3,
                                                                       13920, 3876, 13956, 448,
                                                                       458, 4542, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15072, 0, 3,
                                                                       13992, 3972, 14052, 478,
                                                                       493, 4662, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15162, 0, 3,
                                                                       14052, 4002, 14112, 493,
                                                                       508, 4707, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15252, 0, 3,
                                                                       14112, 4032, 14172, 508,
                                                                       523, 4752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15342, 0, 3,
                                                                       14172, 4062, 14232, 523,
                                                                       538, 4797, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15432, 0, 3,
                                                                       14232, 4092, 14292, 538,
                                                                       553, 4842, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15522, 0, 3,
                                                                       14292, 4122, 14352, 553,
                                                                       568, 4887, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15612, 0, 3,
                                                                       14352, 4152, 14412, 568,
                                                                       583, 4932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15702, 0, 3,
                                                                       14412, 4182, 14472, 583,
                                                                       598, 4977, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15792, 0, 3,
                                                                       14532, 4302, 14592, 628,
                                                                       643, 5112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15882, 0, 3,
                                                                       14592, 4332, 14652, 643,
                                                                       658, 5157, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15972, 0, 3,
                                                                       14652, 4362, 14712, 658,
                                                                       673, 5202, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16062, 0, 3,
                                                                       14712, 4392, 14772, 673,
                                                                       688, 5247, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16152, 0, 3,
                                                                       14772, 4422, 14832, 688,
                                                                       703, 5292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16242, 0, 3,
                                                                       14832, 4452, 14892, 703,
                                                                       718, 5337, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16332, 0, 3,
                                                                       14892, 4482, 14952, 718,
                                                                       733, 5382, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16422, 0, 3,
                                                                       14952, 4512, 15012, 733,
                                                                       748, 5427, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16512, 0, 3,
                                                                       15072, 4662, 15162, 778,
                                                                       799, 5598, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16638, 0, 3,
                                                                       15162, 4707, 15252, 799,
                                                                       820, 5661, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16764, 0, 3,
                                                                       15252, 4752, 15342, 820,
                                                                       841, 5724, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16890, 0, 3,
                                                                       15342, 4797, 15432, 841,
                                                                       862, 5787, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17016, 0, 3,
                                                                       15432, 4842, 15522, 862,
                                                                       883, 5850, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17142, 0, 3,
                                                                       15522, 4887, 15612, 883,
                                                                       904, 5913, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17268, 0, 3,
                                                                       15612, 4932, 15702, 904,
                                                                       925, 5976, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17394, 0, 3,
                                                                       15792, 5112, 15882, 967,
                                                                       988, 6165, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17520, 0, 3,
                                                                       15882, 5157, 15972, 988,
                                                                       1009, 6228, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17646, 0, 3,
                                                                       15972, 5202, 16062, 1009,
                                                                       1030, 6291, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17772, 0, 3,
                                                                       16062, 5247, 16152, 1030,
                                                                       1051, 6354, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17898, 0, 3,
                                                                       16152, 5292, 16242, 1051,
                                                                       1072, 6417, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18024, 0, 3,
                                                                       16242, 5337, 16332, 1072,
                                                                       1093, 6480, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18150, 0, 3,
                                                                       16332, 5382, 16422, 1093,
                                                                       1114, 6543, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18276, 0, 3,
                                                                       16512, 5598, 16638, 1156,
                                                                       1184, 6774, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18444, 0, 3,
                                                                       16638, 5661, 16764, 1184,
                                                                       1212, 6858, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18612, 0, 3,
                                                                       16764, 5724, 16890, 1212,
                                                                       1240, 6942, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18780, 0, 3,
                                                                       16890, 5787, 17016, 1240,
                                                                       1268, 7026, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18948, 0, 3,
                                                                       17016, 5850, 17142, 1268,
                                                                       1296, 7110, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19116, 0, 3,
                                                                       17142, 5913, 17268, 1296,
                                                                       1324, 7194, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19284, 0, 3,
                                                                       17394, 6165, 17520, 1380,
                                                                       1408, 7446, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19452, 0, 3,
                                                                       17520, 6228, 17646, 1408,
                                                                       1436, 7530, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19620, 0, 3,
                                                                       17646, 6291, 17772, 1436,
                                                                       1464, 7614, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19788, 0, 3,
                                                                       17772, 6354, 17898, 1464,
                                                                       1492, 7698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19956, 0, 3,
                                                                       17898, 6417, 18024, 1492,
                                                                       1520, 7782, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20124, 0, 3,
                                                                       18024, 6480, 18150, 1520,
                                                                       1548, 7866, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20292, 0, 3,
                                                                       18276, 6774, 18444, 1604,
                                                                       1640, 8166, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20508, 0, 3,
                                                                       18444, 6858, 18612, 1640,
                                                                       1676, 8274, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20724, 0, 3,
                                                                       18612, 6942, 18780, 1676,
                                                                       1712, 8382, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20940, 0, 3,
                                                                       18780, 7026, 18948, 1712,
                                                                       1748, 8490, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21156, 0, 3,
                                                                       18948, 7110, 19116, 1748,
                                                                       1784, 8598, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21372, 0, 3,
                                                                       19284, 7446, 19452, 1856,
                                                                       1892, 8922, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21588, 0, 3,
                                                                       19452, 7530, 19620, 1892,
                                                                       1928, 9030, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21804, 0, 3,
                                                                       19620, 7614, 19788, 1928,
                                                                       1964, 9138, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22020, 0, 3,
                                                                       19788, 7698, 19956, 1964,
                                                                       2000, 9246, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22236, 0, 3,
                                                                       19956, 7782, 20124, 2000,
                                                                       2036, 9354, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 22452, 0, 3,
                                                                       20292, 8166, 20508, 2108,
                                                                       2153, 9732, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 22722, 0, 3,
                                                                       20508, 8274, 20724, 2153,
                                                                       2198, 9867, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 22992, 0, 3,
                                                                       20724, 8382, 20940, 2198,
                                                                       2243, 10002, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23262, 0, 3,
                                                                       20940, 8490, 21156, 2243,
                                                                       2288, 10137, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23532, 0, 3,
                                                                       21372, 8922, 21588, 2378,
                                                                       2423, 10542, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23802, 0, 3,
                                                                       21588, 9030, 21804, 2423,
                                                                       2468, 10677, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24072, 0, 3,
                                                                       21804, 9138, 22020, 2468,
                                                                       2513, 10812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24342, 0, 3,
                                                                       22020, 9246, 22236, 2513,
                                                                       2558, 10947, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 24612, 0, 3,
                                                                       22452, 9732, 22722, 2648,
                                                                       2703, 11412, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 24942, 0, 3,
                                                                       22722, 9867, 22992, 2703,
                                                                       2758, 11577, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 25272, 0, 3,
                                                                       22992, 10002, 23262, 2758,
                                                                       2813, 11742, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 25602, 0, 3,
                                                                       23532, 10542, 23802, 2923,
                                                                       2978, 12237, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 25932, 0, 3,
                                                                       23802, 10677, 24072, 2978,
                                                                       3033, 12402, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 26262, 0, 3,
                                                                       24072, 10812, 24342, 3033,
                                                                       3088, 12567, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26592, 3, 3198,
                                                                       3201, 12732, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26602, 3, 3201,
                                                                       3204, 12738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26612, 3, 3204,
                                                                       3207, 12744, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26622, 3, 3207,
                                                                       3210, 12750, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26632, 3, 3210,
                                                                       3213, 12756, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26642, 3, 3213,
                                                                       3216, 12762, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26652, 3, 3216,
                                                                       3219, 12768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26662, 3, 3219,
                                                                       3222, 12774, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26672, 3, 3222,
                                                                       3225, 12780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26682, 3, 3225,
                                                                       3228, 12786, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26692, 3, 3228,
                                                                       3231, 12792, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26702, 3, 3231,
                                                                       3234, 12798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26712, 3, 3240,
                                                                       3243, 12804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26722, 3, 3243,
                                                                       3246, 12810, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26732, 3, 3246,
                                                                       3249, 12816, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26742, 3, 3249,
                                                                       3252, 12822, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26752, 3, 3252,
                                                                       3255, 12828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26762, 3, 3255,
                                                                       3258, 12834, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26772, 3, 3258,
                                                                       3261, 12840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26782, 3, 3261,
                                                                       3264, 12846, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26792, 3, 3264,
                                                                       3267, 12852, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26802, 3, 3267,
                                                                       3270, 12858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26812, 3, 3270,
                                                                       3273, 12864, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 26822, 3, 3273,
                                                                       3276, 12870, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26832, 0, 3,
                                                                       26592, 12732, 26602,
                                                                       12876, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26862, 0, 3,
                                                                       26602, 12738, 26612,
                                                                       12894, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26892, 0, 3,
                                                                       26612, 12744, 26622,
                                                                       12912, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26922, 0, 3,
                                                                       26622, 12750, 26632,
                                                                       12930, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26952, 0, 3,
                                                                       26632, 12756, 26642,
                                                                       12948, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 26982, 0, 3,
                                                                       26642, 12762, 26652,
                                                                       12966, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27012, 0, 3,
                                                                       26652, 12768, 26662,
                                                                       12984, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27042, 0, 3,
                                                                       26662, 12774, 26672,
                                                                       13002, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27072, 0, 3,
                                                                       26672, 12780, 26682,
                                                                       13020, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27102, 0, 3,
                                                                       26682, 12786, 26692,
                                                                       13038, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27132, 0, 3,
                                                                       26692, 12792, 26702,
                                                                       13056, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27162, 0, 3,
                                                                       26712, 12804, 26722,
                                                                       13074, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27192, 0, 3,
                                                                       26722, 12810, 26732,
                                                                       13092, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27222, 0, 3,
                                                                       26732, 12816, 26742,
                                                                       13110, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27252, 0, 3,
                                                                       26742, 12822, 26752,
                                                                       13128, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27282, 0, 3,
                                                                       26752, 12828, 26762,
                                                                       13146, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27312, 0, 3,
                                                                       26762, 12834, 26772,
                                                                       13164, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27342, 0, 3,
                                                                       26772, 12840, 26782,
                                                                       13182, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27372, 0, 3,
                                                                       26782, 12846, 26792,
                                                                       13200, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27402, 0, 3,
                                                                       26792, 12852, 26802,
                                                                       13218, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27432, 0, 3,
                                                                       26802, 12858, 26812,
                                                                       13236, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 27462, 0, 3,
                                                                       26812, 12864, 26822,
                                                                       13254, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27492, 0, 3,
                                                                       26832, 12876, 26862, 3480,
                                                                       3498, 13272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27552, 0, 3,
                                                                       26862, 12894, 26892, 3498,
                                                                       3516, 13308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27612, 0, 3,
                                                                       26892, 12912, 26922, 3516,
                                                                       3534, 13344, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27672, 0, 3,
                                                                       26922, 12930, 26952, 3534,
                                                                       3552, 13380, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27732, 0, 3,
                                                                       26952, 12948, 26982, 3552,
                                                                       3570, 13416, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27792, 0, 3,
                                                                       26982, 12966, 27012, 3570,
                                                                       3588, 13452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27852, 0, 3,
                                                                       27012, 12984, 27042, 3588,
                                                                       3606, 13488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27912, 0, 3,
                                                                       27042, 13002, 27072, 3606,
                                                                       3624, 13524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 27972, 0, 3,
                                                                       27072, 13020, 27102, 3624,
                                                                       3642, 13560, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28032, 0, 3,
                                                                       27102, 13038, 27132, 3642,
                                                                       3660, 13596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28092, 0, 3,
                                                                       27162, 13074, 27192, 3696,
                                                                       3714, 13632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28152, 0, 3,
                                                                       27192, 13092, 27222, 3714,
                                                                       3732, 13668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28212, 0, 3,
                                                                       27222, 13110, 27252, 3732,
                                                                       3750, 13704, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28272, 0, 3,
                                                                       27252, 13128, 27282, 3750,
                                                                       3768, 13740, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28332, 0, 3,
                                                                       27282, 13146, 27312, 3768,
                                                                       3786, 13776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28392, 0, 3,
                                                                       27312, 13164, 27342, 3786,
                                                                       3804, 13812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28452, 0, 3,
                                                                       27342, 13182, 27372, 3804,
                                                                       3822, 13848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28512, 0, 3,
                                                                       27372, 13200, 27402, 3822,
                                                                       3840, 13884, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28572, 0, 3,
                                                                       27402, 13218, 27432, 3840,
                                                                       3858, 13920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 28632, 0, 3,
                                                                       27432, 13236, 27462, 3858,
                                                                       3876, 13956, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28692, 0, 3,
                                                                       27492, 13272, 27552, 3912,
                                                                       3942, 13992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28792, 0, 3,
                                                                       27552, 13308, 27612, 3942,
                                                                       3972, 14052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28892, 0, 3,
                                                                       27612, 13344, 27672, 3972,
                                                                       4002, 14112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 28992, 0, 3,
                                                                       27672, 13380, 27732, 4002,
                                                                       4032, 14172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29092, 0, 3,
                                                                       27732, 13416, 27792, 4032,
                                                                       4062, 14232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29192, 0, 3,
                                                                       27792, 13452, 27852, 4062,
                                                                       4092, 14292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29292, 0, 3,
                                                                       27852, 13488, 27912, 4092,
                                                                       4122, 14352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29392, 0, 3,
                                                                       27912, 13524, 27972, 4122,
                                                                       4152, 14412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29492, 0, 3,
                                                                       27972, 13560, 28032, 4152,
                                                                       4182, 14472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29592, 0, 3,
                                                                       28092, 13632, 28152, 4242,
                                                                       4272, 14532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29692, 0, 3,
                                                                       28152, 13668, 28212, 4272,
                                                                       4302, 14592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29792, 0, 3,
                                                                       28212, 13704, 28272, 4302,
                                                                       4332, 14652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29892, 0, 3,
                                                                       28272, 13740, 28332, 4332,
                                                                       4362, 14712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 29992, 0, 3,
                                                                       28332, 13776, 28392, 4362,
                                                                       4392, 14772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30092, 0, 3,
                                                                       28392, 13812, 28452, 4392,
                                                                       4422, 14832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30192, 0, 3,
                                                                       28452, 13848, 28512, 4422,
                                                                       4452, 14892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30292, 0, 3,
                                                                       28512, 13884, 28572, 4452,
                                                                       4482, 14952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30392, 0, 3,
                                                                       28572, 13920, 28632, 4482,
                                                                       4512, 15012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30492, 0, 3,
                                                                       28692, 13992, 28792, 4572,
                                                                       4617, 15072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30642, 0, 3,
                                                                       28792, 14052, 28892, 4617,
                                                                       4662, 15162, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30792, 0, 3,
                                                                       28892, 14112, 28992, 4662,
                                                                       4707, 15252, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30942, 0, 3,
                                                                       28992, 14172, 29092, 4707,
                                                                       4752, 15342, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31092, 0, 3,
                                                                       29092, 14232, 29192, 4752,
                                                                       4797, 15432, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31242, 0, 3,
                                                                       29192, 14292, 29292, 4797,
                                                                       4842, 15522, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31392, 0, 3,
                                                                       29292, 14352, 29392, 4842,
                                                                       4887, 15612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31542, 0, 3,
                                                                       29392, 14412, 29492, 4887,
                                                                       4932, 15702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31692, 0, 3,
                                                                       29592, 14532, 29692, 5022,
                                                                       5067, 15792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31842, 0, 3,
                                                                       29692, 14592, 29792, 5067,
                                                                       5112, 15882, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31992, 0, 3,
                                                                       29792, 14652, 29892, 5112,
                                                                       5157, 15972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32142, 0, 3,
                                                                       29892, 14712, 29992, 5157,
                                                                       5202, 16062, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32292, 0, 3,
                                                                       29992, 14772, 30092, 5202,
                                                                       5247, 16152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32442, 0, 3,
                                                                       30092, 14832, 30192, 5247,
                                                                       5292, 16242, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32592, 0, 3,
                                                                       30192, 14892, 30292, 5292,
                                                                       5337, 16332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32742, 0, 3,
                                                                       30292, 14952, 30392, 5337,
                                                                       5382, 16422, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 32892, 0, 3,
                                                                       30492, 15072, 30642, 5472,
                                                                       5535, 16512, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33102, 0, 3,
                                                                       30642, 15162, 30792, 5535,
                                                                       5598, 16638, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33312, 0, 3,
                                                                       30792, 15252, 30942, 5598,
                                                                       5661, 16764, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33522, 0, 3,
                                                                       30942, 15342, 31092, 5661,
                                                                       5724, 16890, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33732, 0, 3,
                                                                       31092, 15432, 31242, 5724,
                                                                       5787, 17016, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33942, 0, 3,
                                                                       31242, 15522, 31392, 5787,
                                                                       5850, 17142, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34152, 0, 3,
                                                                       31392, 15612, 31542, 5850,
                                                                       5913, 17268, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34362, 0, 3,
                                                                       31692, 15792, 31842, 6039,
                                                                       6102, 17394, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34572, 0, 3,
                                                                       31842, 15882, 31992, 6102,
                                                                       6165, 17520, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34782, 0, 3,
                                                                       31992, 15972, 32142, 6165,
                                                                       6228, 17646, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34992, 0, 3,
                                                                       32142, 16062, 32292, 6228,
                                                                       6291, 17772, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35202, 0, 3,
                                                                       32292, 16152, 32442, 6291,
                                                                       6354, 17898, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35412, 0, 3,
                                                                       32442, 16242, 32592, 6354,
                                                                       6417, 18024, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35622, 0, 3,
                                                                       32592, 16332, 32742, 6417,
                                                                       6480, 18150, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 35832, 0, 3,
                                                                       32892, 16512, 33102, 6606,
                                                                       6690, 18276, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36112, 0, 3,
                                                                       33102, 16638, 33312, 6690,
                                                                       6774, 18444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36392, 0, 3,
                                                                       33312, 16764, 33522, 6774,
                                                                       6858, 18612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36672, 0, 3,
                                                                       33522, 16890, 33732, 6858,
                                                                       6942, 18780, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36952, 0, 3,
                                                                       33732, 17016, 33942, 6942,
                                                                       7026, 18948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 37232, 0, 3,
                                                                       33942, 17142, 34152, 7026,
                                                                       7110, 19116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 37512, 0, 3,
                                                                       34362, 17394, 34572, 7278,
                                                                       7362, 19284, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 37792, 0, 3,
                                                                       34572, 17520, 34782, 7362,
                                                                       7446, 19452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38072, 0, 3,
                                                                       34782, 17646, 34992, 7446,
                                                                       7530, 19620, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38352, 0, 3,
                                                                       34992, 17772, 35202, 7530,
                                                                       7614, 19788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38632, 0, 3,
                                                                       35202, 17898, 35412, 7614,
                                                                       7698, 19956, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38912, 0, 3,
                                                                       35412, 18024, 35622, 7698,
                                                                       7782, 20124, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 39192, 0, 3,
                                                                       35832, 18276, 36112, 7950,
                                                                       8058, 20292, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 39552, 0, 3,
                                                                       36112, 18444, 36392, 8058,
                                                                       8166, 20508, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 39912, 0, 3,
                                                                       36392, 18612, 36672, 8166,
                                                                       8274, 20724, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 40272, 0, 3,
                                                                       36672, 18780, 36952, 8274,
                                                                       8382, 20940, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 40632, 0, 3,
                                                                       36952, 18948, 37232, 8382,
                                                                       8490, 21156, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 40992, 0, 3,
                                                                       37512, 19284, 37792, 8706,
                                                                       8814, 21372, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41352, 0, 3,
                                                                       37792, 19452, 38072, 8814,
                                                                       8922, 21588, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41712, 0, 3,
                                                                       38072, 19620, 38352, 8922,
                                                                       9030, 21804, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42072, 0, 3,
                                                                       38352, 19788, 38632, 9030,
                                                                       9138, 22020, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42432, 0, 3,
                                                                       38632, 19956, 38912, 9138,
                                                                       9246, 22236, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 42792, 0, 3,
                                                                       39192, 20292, 39552, 9462,
                                                                       9597, 22452, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 43242, 0, 3,
                                                                       39552, 20508, 39912, 9597,
                                                                       9732, 22722, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 43692, 0, 3,
                                                                       39912, 20724, 40272, 9732,
                                                                       9867, 22992, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 44142, 0, 3,
                                                                       40272, 20940, 40632, 9867,
                                                                       10002, 23262, ncols,
                                                                       gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 44592, 0, 3,
                                                                       40992, 21372, 41352,
                                                                       10272, 10407, 23532,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 45042, 0, 3,
                                                                       41352, 21588, 41712,
                                                                       10407, 10542, 23802,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 45492, 0, 3,
                                                                       41712, 21804, 42072,
                                                                       10542, 10677, 24072,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 45942, 0, 3,
                                                                       42072, 22020, 42432,
                                                                       10677, 10812, 24342,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 46392, 0, 3,
                                                                       42792, 22452, 43242,
                                                                       11082, 11247, 24612,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 46942, 0, 3,
                                                                       43242, 22722, 43692,
                                                                       11247, 11412, 24942,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 47492, 0, 3,
                                                                       43692, 22992, 44142,
                                                                       11412, 11577, 25272,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 48042, 0, 3,
                                                                       44592, 23532, 45042,
                                                                       11907, 12072, 25602,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 48592, 0, 3,
                                                                       45042, 23802, 45492,
                                                                       12072, 12237, 25932,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 49142, 0, 3,
                                                                       45492, 24072, 45942,
                                                                       12237, 12402, 26262,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49692, 3, 12732,
                                                                       12738, 26612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49707, 3, 12738,
                                                                       12744, 26622, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49722, 3, 12744,
                                                                       12750, 26632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49737, 3, 12750,
                                                                       12756, 26642, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49752, 3, 12756,
                                                                       12762, 26652, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49767, 3, 12762,
                                                                       12768, 26662, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49782, 3, 12768,
                                                                       12774, 26672, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49797, 3, 12774,
                                                                       12780, 26682, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49812, 3, 12780,
                                                                       12786, 26692, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49827, 3, 12786,
                                                                       12792, 26702, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49842, 3, 12804,
                                                                       12810, 26732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49857, 3, 12810,
                                                                       12816, 26742, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49872, 3, 12816,
                                                                       12822, 26752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49887, 3, 12822,
                                                                       12828, 26762, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49902, 3, 12828,
                                                                       12834, 26772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49917, 3, 12834,
                                                                       12840, 26782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49932, 3, 12840,
                                                                       12846, 26792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49947, 3, 12846,
                                                                       12852, 26802, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49962, 3, 12852,
                                                                       12858, 26812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49977, 3, 12858,
                                                                       12864, 26822, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 49992, 0, 3,
                                                                       49692, 26612, 49707,
                                                                       12876, 12894, 26892,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50037, 0, 3,
                                                                       49707, 26622, 49722,
                                                                       12894, 12912, 26922,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50082, 0, 3,
                                                                       49722, 26632, 49737,
                                                                       12912, 12930, 26952,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50127, 0, 3,
                                                                       49737, 26642, 49752,
                                                                       12930, 12948, 26982,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50172, 0, 3,
                                                                       49752, 26652, 49767,
                                                                       12948, 12966, 27012,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50217, 0, 3,
                                                                       49767, 26662, 49782,
                                                                       12966, 12984, 27042,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50262, 0, 3,
                                                                       49782, 26672, 49797,
                                                                       12984, 13002, 27072,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50307, 0, 3,
                                                                       49797, 26682, 49812,
                                                                       13002, 13020, 27102,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50352, 0, 3,
                                                                       49812, 26692, 49827,
                                                                       13020, 13038, 27132,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50397, 0, 3,
                                                                       49842, 26732, 49857,
                                                                       13074, 13092, 27222,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50442, 0, 3,
                                                                       49857, 26742, 49872,
                                                                       13092, 13110, 27252,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50487, 0, 3,
                                                                       49872, 26752, 49887,
                                                                       13110, 13128, 27282,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50532, 0, 3,
                                                                       49887, 26762, 49902,
                                                                       13128, 13146, 27312,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50577, 0, 3,
                                                                       49902, 26772, 49917,
                                                                       13146, 13164, 27342,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50622, 0, 3,
                                                                       49917, 26782, 49932,
                                                                       13164, 13182, 27372,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50667, 0, 3,
                                                                       49932, 26792, 49947,
                                                                       13182, 13200, 27402,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50712, 0, 3,
                                                                       49947, 26802, 49962,
                                                                       13200, 13218, 27432,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50757, 0, 3,
                                                                       49962, 26812, 49977,
                                                                       13218, 13236, 27462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50802, 0, 3,
                                                                       49992, 26892, 50037,
                                                                       13272, 13308, 27612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50892, 0, 3,
                                                                       50037, 26922, 50082,
                                                                       13308, 13344, 27672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50982, 0, 3,
                                                                       50082, 26952, 50127,
                                                                       13344, 13380, 27732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51072, 0, 3,
                                                                       50127, 26982, 50172,
                                                                       13380, 13416, 27792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51162, 0, 3,
                                                                       50172, 27012, 50217,
                                                                       13416, 13452, 27852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51252, 0, 3,
                                                                       50217, 27042, 50262,
                                                                       13452, 13488, 27912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51342, 0, 3,
                                                                       50262, 27072, 50307,
                                                                       13488, 13524, 27972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51432, 0, 3,
                                                                       50307, 27102, 50352,
                                                                       13524, 13560, 28032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51522, 0, 3,
                                                                       50397, 27222, 50442,
                                                                       13632, 13668, 28212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51612, 0, 3,
                                                                       50442, 27252, 50487,
                                                                       13668, 13704, 28272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51702, 0, 3,
                                                                       50487, 27282, 50532,
                                                                       13704, 13740, 28332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51792, 0, 3,
                                                                       50532, 27312, 50577,
                                                                       13740, 13776, 28392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51882, 0, 3,
                                                                       50577, 27342, 50622,
                                                                       13776, 13812, 28452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51972, 0, 3,
                                                                       50622, 27372, 50667,
                                                                       13812, 13848, 28512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 52062, 0, 3,
                                                                       50667, 27402, 50712,
                                                                       13848, 13884, 28572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 52152, 0, 3,
                                                                       50712, 27432, 50757,
                                                                       13884, 13920, 28632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52242, 0, 3,
                                                                       50802, 27612, 50892,
                                                                       13992, 14052, 28892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52392, 0, 3,
                                                                       50892, 27672, 50982,
                                                                       14052, 14112, 28992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52542, 0, 3,
                                                                       50982, 27732, 51072,
                                                                       14112, 14172, 29092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52692, 0, 3,
                                                                       51072, 27792, 51162,
                                                                       14172, 14232, 29192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52842, 0, 3,
                                                                       51162, 27852, 51252,
                                                                       14232, 14292, 29292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52992, 0, 3,
                                                                       51252, 27912, 51342,
                                                                       14292, 14352, 29392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53142, 0, 3,
                                                                       51342, 27972, 51432,
                                                                       14352, 14412, 29492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53292, 0, 3,
                                                                       51522, 28212, 51612,
                                                                       14532, 14592, 29792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53442, 0, 3,
                                                                       51612, 28272, 51702,
                                                                       14592, 14652, 29892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53592, 0, 3,
                                                                       51702, 28332, 51792,
                                                                       14652, 14712, 29992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53742, 0, 3,
                                                                       51792, 28392, 51882,
                                                                       14712, 14772, 30092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53892, 0, 3,
                                                                       51882, 28452, 51972,
                                                                       14772, 14832, 30192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 54042, 0, 3,
                                                                       51972, 28512, 52062,
                                                                       14832, 14892, 30292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 54192, 0, 3,
                                                                       52062, 28572, 52152,
                                                                       14892, 14952, 30392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54342, 0, 3,
                                                                       52242, 28892, 52392,
                                                                       15072, 15162, 30792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54567, 0, 3,
                                                                       52392, 28992, 52542,
                                                                       15162, 15252, 30942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54792, 0, 3,
                                                                       52542, 29092, 52692,
                                                                       15252, 15342, 31092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55017, 0, 3,
                                                                       52692, 29192, 52842,
                                                                       15342, 15432, 31242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55242, 0, 3,
                                                                       52842, 29292, 52992,
                                                                       15432, 15522, 31392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55467, 0, 3,
                                                                       52992, 29392, 53142,
                                                                       15522, 15612, 31542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55692, 0, 3,
                                                                       53292, 29792, 53442,
                                                                       15792, 15882, 31992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55917, 0, 3,
                                                                       53442, 29892, 53592,
                                                                       15882, 15972, 32142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 56142, 0, 3,
                                                                       53592, 29992, 53742,
                                                                       15972, 16062, 32292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 56367, 0, 3,
                                                                       53742, 30092, 53892,
                                                                       16062, 16152, 32442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 56592, 0, 3,
                                                                       53892, 30192, 54042,
                                                                       16152, 16242, 32592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 56817, 0, 3,
                                                                       54042, 30292, 54192,
                                                                       16242, 16332, 32742,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 57042, 0, 3,
                                                                       54342, 30792, 54567,
                                                                       16512, 16638, 33312,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 57357, 0, 3,
                                                                       54567, 30942, 54792,
                                                                       16638, 16764, 33522,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 57672, 0, 3,
                                                                       54792, 31092, 55017,
                                                                       16764, 16890, 33732,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 57987, 0, 3,
                                                                       55017, 31242, 55242,
                                                                       16890, 17016, 33942,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 58302, 0, 3,
                                                                       55242, 31392, 55467,
                                                                       17016, 17142, 34152,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 58617, 0, 3,
                                                                       55692, 31992, 55917,
                                                                       17394, 17520, 34782,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 58932, 0, 3,
                                                                       55917, 32142, 56142,
                                                                       17520, 17646, 34992,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 59247, 0, 3,
                                                                       56142, 32292, 56367,
                                                                       17646, 17772, 35202,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 59562, 0, 3,
                                                                       56367, 32442, 56592,
                                                                       17772, 17898, 35412,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 59877, 0, 3,
                                                                       56592, 32592, 56817,
                                                                       17898, 18024, 35622,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 60192, 0, 3,
                                                                       57042, 33312, 57357,
                                                                       18276, 18444, 36392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 60612, 0, 3,
                                                                       57357, 33522, 57672,
                                                                       18444, 18612, 36672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 61032, 0, 3,
                                                                       57672, 33732, 57987,
                                                                       18612, 18780, 36952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 61452, 0, 3,
                                                                       57987, 33942, 58302,
                                                                       18780, 18948, 37232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 61872, 0, 3,
                                                                       58617, 34782, 58932,
                                                                       19284, 19452, 38072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 62292, 0, 3,
                                                                       58932, 34992, 59247,
                                                                       19452, 19620, 38352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 62712, 0, 3,
                                                                       59247, 35202, 59562,
                                                                       19620, 19788, 38632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 63132, 0, 3,
                                                                       59562, 35412, 59877,
                                                                       19788, 19956, 38912,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 63552, 0, 3,
                                                                       60192, 36392, 60612,
                                                                       20292, 20508, 39912,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 64092, 0, 3,
                                                                       60612, 36672, 61032,
                                                                       20508, 20724, 40272,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 64632, 0, 3,
                                                                       61032, 36952, 61452,
                                                                       20724, 20940, 40632,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 65172, 0, 3,
                                                                       61872, 38072, 62292,
                                                                       21372, 21588, 41712,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 65712, 0, 3,
                                                                       62292, 38352, 62712,
                                                                       21588, 21804, 42072,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 66252, 0, 3,
                                                                       62712, 38632, 63132,
                                                                       21804, 22020, 42432,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 66792, 0, 3,
                                                                       63552, 39912, 64092,
                                                                       22452, 22722, 43692,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 67467, 0, 3,
                                                                       64092, 40272, 64632,
                                                                       22722, 22992, 44142,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 68142, 0, 3,
                                                                       65172, 41712, 65712,
                                                                       23532, 23802, 45492,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 68817, 0, 3,
                                                                       65712, 42072, 66252,
                                                                       23802, 24072, 45942,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 69492, 0, 3,
                                                                       66792, 43692, 67467,
                                                                       24612, 24942, 47492,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 70317, 0, 3,
                                                                       68142, 45492, 68817,
                                                                       25602, 25932, 49142,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71142, 3, 26592,
                                                                       26602, 49692, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71163, 3, 26602,
                                                                       26612, 49707, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71184, 3, 26612,
                                                                       26622, 49722, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71205, 3, 26622,
                                                                       26632, 49737, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71226, 3, 26632,
                                                                       26642, 49752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71247, 3, 26642,
                                                                       26652, 49767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71268, 3, 26652,
                                                                       26662, 49782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71289, 3, 26662,
                                                                       26672, 49797, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71310, 3, 26672,
                                                                       26682, 49812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71331, 3, 26682,
                                                                       26692, 49827, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71352, 3, 26712,
                                                                       26722, 49842, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71373, 3, 26722,
                                                                       26732, 49857, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71394, 3, 26732,
                                                                       26742, 49872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71415, 3, 26742,
                                                                       26752, 49887, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71436, 3, 26752,
                                                                       26762, 49902, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71457, 3, 26762,
                                                                       26772, 49917, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71478, 3, 26772,
                                                                       26782, 49932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71499, 3, 26782,
                                                                       26792, 49947, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71520, 3, 26792,
                                                                       26802, 49962, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 71541, 3, 26802,
                                                                       26812, 49977, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 71562, 0, 3,
                                                                       71142, 49692, 71163,
                                                                       26832, 26862, 49992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 71625, 0, 3,
                                                                       71163, 49707, 71184,
                                                                       26862, 26892, 50037,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 71688, 0, 3,
                                                                       71184, 49722, 71205,
                                                                       26892, 26922, 50082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 71751, 0, 3,
                                                                       71205, 49737, 71226,
                                                                       26922, 26952, 50127,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 71814, 0, 3,
                                                                       71226, 49752, 71247,
                                                                       26952, 26982, 50172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 71877, 0, 3,
                                                                       71247, 49767, 71268,
                                                                       26982, 27012, 50217,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 71940, 0, 3,
                                                                       71268, 49782, 71289,
                                                                       27012, 27042, 50262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72003, 0, 3,
                                                                       71289, 49797, 71310,
                                                                       27042, 27072, 50307,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72066, 0, 3,
                                                                       71310, 49812, 71331,
                                                                       27072, 27102, 50352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72129, 0, 3,
                                                                       71352, 49842, 71373,
                                                                       27162, 27192, 50397,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72192, 0, 3,
                                                                       71373, 49857, 71394,
                                                                       27192, 27222, 50442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72255, 0, 3,
                                                                       71394, 49872, 71415,
                                                                       27222, 27252, 50487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72318, 0, 3,
                                                                       71415, 49887, 71436,
                                                                       27252, 27282, 50532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72381, 0, 3,
                                                                       71436, 49902, 71457,
                                                                       27282, 27312, 50577,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72444, 0, 3,
                                                                       71457, 49917, 71478,
                                                                       27312, 27342, 50622,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72507, 0, 3,
                                                                       71478, 49932, 71499,
                                                                       27342, 27372, 50667,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72570, 0, 3,
                                                                       71499, 49947, 71520,
                                                                       27372, 27402, 50712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 72633, 0, 3,
                                                                       71520, 49962, 71541,
                                                                       27402, 27432, 50757,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 72696, 0, 3,
                                                                       71562, 49992, 71625,
                                                                       27492, 27552, 50802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 72822, 0, 3,
                                                                       71625, 50037, 71688,
                                                                       27552, 27612, 50892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 72948, 0, 3,
                                                                       71688, 50082, 71751,
                                                                       27612, 27672, 50982,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 73074, 0, 3,
                                                                       71751, 50127, 71814,
                                                                       27672, 27732, 51072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 73200, 0, 3,
                                                                       71814, 50172, 71877,
                                                                       27732, 27792, 51162,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 73326, 0, 3,
                                                                       71877, 50217, 71940,
                                                                       27792, 27852, 51252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 73452, 0, 3,
                                                                       71940, 50262, 72003,
                                                                       27852, 27912, 51342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 73578, 0, 3,
                                                                       72003, 50307, 72066,
                                                                       27912, 27972, 51432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 73704, 0, 3,
                                                                       72129, 50397, 72192,
                                                                       28092, 28152, 51522,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 73830, 0, 3,
                                                                       72192, 50442, 72255,
                                                                       28152, 28212, 51612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 73956, 0, 3,
                                                                       72255, 50487, 72318,
                                                                       28212, 28272, 51702,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 74082, 0, 3,
                                                                       72318, 50532, 72381,
                                                                       28272, 28332, 51792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 74208, 0, 3,
                                                                       72381, 50577, 72444,
                                                                       28332, 28392, 51882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 74334, 0, 3,
                                                                       72444, 50622, 72507,
                                                                       28392, 28452, 51972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 74460, 0, 3,
                                                                       72507, 50667, 72570,
                                                                       28452, 28512, 52062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 74586, 0, 3,
                                                                       72570, 50712, 72633,
                                                                       28512, 28572, 52152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 74712, 0, 3,
                                                                       72696, 50802, 72822,
                                                                       28692, 28792, 52242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 74922, 0, 3,
                                                                       72822, 50892, 72948,
                                                                       28792, 28892, 52392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 75132, 0, 3,
                                                                       72948, 50982, 73074,
                                                                       28892, 28992, 52542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 75342, 0, 3,
                                                                       73074, 51072, 73200,
                                                                       28992, 29092, 52692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 75552, 0, 3,
                                                                       73200, 51162, 73326,
                                                                       29092, 29192, 52842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 75762, 0, 3,
                                                                       73326, 51252, 73452,
                                                                       29192, 29292, 52992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 75972, 0, 3,
                                                                       73452, 51342, 73578,
                                                                       29292, 29392, 53142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 76182, 0, 3,
                                                                       73704, 51522, 73830,
                                                                       29592, 29692, 53292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 76392, 0, 3,
                                                                       73830, 51612, 73956,
                                                                       29692, 29792, 53442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 76602, 0, 3,
                                                                       73956, 51702, 74082,
                                                                       29792, 29892, 53592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 76812, 0, 3,
                                                                       74082, 51792, 74208,
                                                                       29892, 29992, 53742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 77022, 0, 3,
                                                                       74208, 51882, 74334,
                                                                       29992, 30092, 53892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 77232, 0, 3,
                                                                       74334, 51972, 74460,
                                                                       30092, 30192, 54042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 77442, 0, 3,
                                                                       74460, 52062, 74586,
                                                                       30192, 30292, 54192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 77652, 0, 3,
                                                                       74712, 52242, 74922,
                                                                       30492, 30642, 54342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 77967, 0, 3,
                                                                       74922, 52392, 75132,
                                                                       30642, 30792, 54567,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 78282, 0, 3,
                                                                       75132, 52542, 75342,
                                                                       30792, 30942, 54792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 78597, 0, 3,
                                                                       75342, 52692, 75552,
                                                                       30942, 31092, 55017,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 78912, 0, 3,
                                                                       75552, 52842, 75762,
                                                                       31092, 31242, 55242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 79227, 0, 3,
                                                                       75762, 52992, 75972,
                                                                       31242, 31392, 55467,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 79542, 0, 3,
                                                                       76182, 53292, 76392,
                                                                       31692, 31842, 55692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 79857, 0, 3,
                                                                       76392, 53442, 76602,
                                                                       31842, 31992, 55917,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 80172, 0, 3,
                                                                       76602, 53592, 76812,
                                                                       31992, 32142, 56142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 80487, 0, 3,
                                                                       76812, 53742, 77022,
                                                                       32142, 32292, 56367,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 80802, 0, 3,
                                                                       77022, 53892, 77232,
                                                                       32292, 32442, 56592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 81117, 0, 3,
                                                                       77232, 54042, 77442,
                                                                       32442, 32592, 56817,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 81432, 0, 3,
                                                                       77652, 54342, 77967,
                                                                       32892, 33102, 57042,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 81873, 0, 3,
                                                                       77967, 54567, 78282,
                                                                       33102, 33312, 57357,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 82314, 0, 3,
                                                                       78282, 54792, 78597,
                                                                       33312, 33522, 57672,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 82755, 0, 3,
                                                                       78597, 55017, 78912,
                                                                       33522, 33732, 57987,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 83196, 0, 3,
                                                                       78912, 55242, 79227,
                                                                       33732, 33942, 58302,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 83637, 0, 3,
                                                                       79542, 55692, 79857,
                                                                       34362, 34572, 58617,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 84078, 0, 3,
                                                                       79857, 55917, 80172,
                                                                       34572, 34782, 58932,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 84519, 0, 3,
                                                                       80172, 56142, 80487,
                                                                       34782, 34992, 59247,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 84960, 0, 3,
                                                                       80487, 56367, 80802,
                                                                       34992, 35202, 59562,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 85401, 0, 3,
                                                                       80802, 56592, 81117,
                                                                       35202, 35412, 59877,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 85842, 0, 3,
                                                                       81432, 57042, 81873,
                                                                       35832, 36112, 60192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 86430, 0, 3,
                                                                       81873, 57357, 82314,
                                                                       36112, 36392, 60612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 87018, 0, 3,
                                                                       82314, 57672, 82755,
                                                                       36392, 36672, 61032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 87606, 0, 3,
                                                                       82755, 57987, 83196,
                                                                       36672, 36952, 61452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 88194, 0, 3,
                                                                       83637, 58617, 84078,
                                                                       37512, 37792, 61872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 88782, 0, 3,
                                                                       84078, 58932, 84519,
                                                                       37792, 38072, 62292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 89370, 0, 3,
                                                                       84519, 59247, 84960,
                                                                       38072, 38352, 62712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 89958, 0, 3,
                                                                       84960, 59562, 85401,
                                                                       38352, 38632, 63132,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 90546, 0, 3,
                                                                       85842, 60192, 86430,
                                                                       39192, 39552, 63552,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 91302, 0, 3,
                                                                       86430, 60612, 87018,
                                                                       39552, 39912, 64092,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 92058, 0, 3,
                                                                       87018, 61032, 87606,
                                                                       39912, 40272, 64632,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 92814, 0, 3,
                                                                       88194, 61872, 88782,
                                                                       40992, 41352, 65172,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 93570, 0, 3,
                                                                       88782, 62292, 89370,
                                                                       41352, 41712, 65712,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 94326, 0, 3,
                                                                       89370, 62712, 89958,
                                                                       41712, 42072, 66252,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 95082, 0, 3,
                                                                       90546, 63552, 91302,
                                                                       42792, 43242, 66792,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 96027, 0, 3,
                                                                       91302, 64092, 92058,
                                                                       43242, 43692, 67467,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 96972, 0, 3,
                                                                       92814, 65172, 93570,
                                                                       44592, 45042, 68142,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 97917, 0, 3,
                                                                       93570, 65712, 94326,
                                                                       45042, 45492, 68817,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 98862, 0, 3,
                                                                       95082, 66792, 96027,
                                                                       46392, 46942, 69492,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 100017, 0, 3,
                                                                       96972, 68142, 97917,
                                                                       48042, 48592, 70317,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 101172, 85842, 588, ncols);

                    simdfunc::contract_primitives(buffer, 102068, 88194, 588, ncols);

                    simdfunc::contract_primitives(buffer, 102964, 90546, 756, ncols);

                    simdfunc::contract_primitives(buffer, 104116, 92814, 756, ncols);

                    simdfunc::contract_primitives(buffer, 105268, 95082, 945, ncols);

                    simdfunc::contract_primitives(buffer, 106708, 96972, 945, ncols);

                    simdfunc::contract_primitives(buffer, 108148, 98862, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 109908, 100017, 1155, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 101760, 101172, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 102656, 102068, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 103720, 102964, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 104872, 104116, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 106213, 105268, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 107653, 106708, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 109303, 108148, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 111063, 109908, 55, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 111668, 101760, 103720, 11, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 112592, 102656, 104872, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 113516, 103720, 106213, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 114704, 104872, 107653, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 115892, 106213, 109303, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 117377, 107653, 111063, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 118862, 111668, 113516, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 120710, 112592, 114704, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 122558, 113516, 115892, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 124934, 114704, 117377, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 127310, 118862, 122558, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 130390, 120710, 124934, 11, nmax);

        simdtrf::transform_i_inner(buffer, 133470, 130390, 10, 11, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 133470, 143, nmax);

        simdtrf::transform_i_inner(buffer, 133470, 127310, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 1001 * nvalues + n * npairs, nvalues, buffer, 133470,
                                   143, nmax);
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
