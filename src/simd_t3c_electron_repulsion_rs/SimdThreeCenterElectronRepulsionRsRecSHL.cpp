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


#include "SimdThreeCenterElectronRepulsionRsRecSHL.hpp"

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
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_shl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_shl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 72955, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 374 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 72955, 70708, 1890, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 13,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 21, 3, 13,
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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1156, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1159, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1162, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1165, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1168, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1171, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1174, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1177, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1180, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1183, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1186, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1189, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1192, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1195, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1198, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1201, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1204, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1207, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1210, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1213, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1216, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1219, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1222, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1225, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1228, 3, 9, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1237, 3, 10, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1246, 3, 11, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1255, 3, 12, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1264, 3, 13, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1273, 3, 14, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1282, 3, 15, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1291, 3, 16, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1300, 3, 17, 66,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1309, 3, 18, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1318, 3, 19, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1327, 3, 24, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1336, 3, 25, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1345, 3, 26, 87,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1354, 3, 27, 90,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1363, 3, 28, 93,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1372, 3, 29, 96,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1381, 3, 30, 99,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1390, 3, 31, 102,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1399, 3, 32, 105,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1408, 3, 33, 108,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1417, 3, 34, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1426, 3, 42, 126,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1444, 3, 45, 132,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1462, 3, 48, 138,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1480, 3, 51, 144,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1498, 3, 54, 150,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1516, 3, 57, 156,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1534, 3, 60, 162,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1552, 3, 63, 168,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1570, 3, 66, 174,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1588, 3, 69, 180,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1606, 3, 81, 198,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1624, 3, 84, 204,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1642, 3, 87, 210,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1660, 3, 90, 216,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1678, 3, 93, 222,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1696, 3, 96, 228,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1714, 3, 99, 234,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1732, 3, 102, 240,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1750, 3, 105, 246,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1768, 3, 108, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1786, 3, 126, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1816, 3, 132, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1846, 3, 138, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1876, 3, 144, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1906, 3, 150, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1936, 3, 156, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1966, 3, 162, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1996, 3, 168, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2026, 3, 174, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2056, 3, 198, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2086, 3, 204, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2116, 3, 210, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2146, 3, 216, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2176, 3, 222, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2206, 3, 228, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2236, 3, 234, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2266, 3, 240, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2296, 3, 246, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2326, 3, 278, 508,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2371, 3, 288, 523,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2416, 3, 298, 538,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2461, 3, 308, 553,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2506, 3, 318, 568,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2551, 3, 328, 583,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2596, 3, 338, 598,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2641, 3, 348, 613,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2686, 3, 388, 658,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2731, 3, 398, 673,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2776, 3, 408, 688,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2821, 3, 418, 703,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2866, 3, 428, 718,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2911, 3, 438, 733,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2956, 3, 448, 748,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3001, 3, 458, 763,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3046, 3, 508, 820,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3109, 3, 523, 841,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3172, 3, 538, 862,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3235, 3, 553, 883,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3298, 3, 568, 904,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3361, 3, 583, 925,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3424, 3, 598, 946,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3487, 3, 658,
                                                                       1009, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3550, 3, 673,
                                                                       1030, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3613, 3, 688,
                                                                       1051, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3676, 3, 703,
                                                                       1072, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3739, 3, 718,
                                                                       1093, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3802, 3, 733,
                                                                       1114, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3865, 3, 748,
                                                                       1135, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3928, 3, 7, 8,
                                                                       1156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3934, 3, 8, 9,
                                                                       1159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3940, 3, 9, 10,
                                                                       1162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3946, 3, 10, 11,
                                                                       1165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3952, 3, 11, 12,
                                                                       1168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3958, 3, 12, 13,
                                                                       1171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3964, 3, 13, 14,
                                                                       1174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3970, 3, 14, 15,
                                                                       1177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3976, 3, 15, 16,
                                                                       1180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3982, 3, 16, 17,
                                                                       1183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3988, 3, 17, 18,
                                                                       1186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3994, 3, 18, 19,
                                                                       1189, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4000, 3, 22, 23,
                                                                       1192, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4006, 3, 23, 24,
                                                                       1195, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4012, 3, 24, 25,
                                                                       1198, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4018, 3, 25, 26,
                                                                       1201, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4024, 3, 26, 27,
                                                                       1204, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4030, 3, 27, 28,
                                                                       1207, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4036, 3, 28, 29,
                                                                       1210, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4042, 3, 29, 30,
                                                                       1213, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4048, 3, 30, 31,
                                                                       1216, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4054, 3, 31, 32,
                                                                       1219, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4060, 3, 32, 33,
                                                                       1222, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4066, 3, 33, 34,
                                                                       1225, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4072, 0, 3, 3928,
                                                                       1156, 3934, 1228, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4090, 0, 3, 3934,
                                                                       1159, 3940, 1237, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4108, 0, 3, 3940,
                                                                       1162, 3946, 1246, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4126, 0, 3, 3946,
                                                                       1165, 3952, 1255, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4144, 0, 3, 3952,
                                                                       1168, 3958, 1264, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4162, 0, 3, 3958,
                                                                       1171, 3964, 1273, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4180, 0, 3, 3964,
                                                                       1174, 3970, 1282, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4198, 0, 3, 3970,
                                                                       1177, 3976, 1291, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4216, 0, 3, 3976,
                                                                       1180, 3982, 1300, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4234, 0, 3, 3982,
                                                                       1183, 3988, 1309, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4252, 0, 3, 3988,
                                                                       1186, 3994, 1318, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4270, 0, 3, 4000,
                                                                       1192, 4006, 1327, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4288, 0, 3, 4006,
                                                                       1195, 4012, 1336, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4306, 0, 3, 4012,
                                                                       1198, 4018, 1345, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4324, 0, 3, 4018,
                                                                       1201, 4024, 1354, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4342, 0, 3, 4024,
                                                                       1204, 4030, 1363, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4360, 0, 3, 4030,
                                                                       1207, 4036, 1372, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4378, 0, 3, 4036,
                                                                       1210, 4042, 1381, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4396, 0, 3, 4042,
                                                                       1213, 4048, 1390, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4414, 0, 3, 4048,
                                                                       1216, 4054, 1399, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4432, 0, 3, 4054,
                                                                       1219, 4060, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4450, 0, 3, 4060,
                                                                       1222, 4066, 1417, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4468, 0, 3, 4072,
                                                                       1228, 4090, 114, 120,
                                                                       1426, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4504, 0, 3, 4090,
                                                                       1237, 4108, 120, 126,
                                                                       1444, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4540, 0, 3, 4108,
                                                                       1246, 4126, 126, 132,
                                                                       1462, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4576, 0, 3, 4126,
                                                                       1255, 4144, 132, 138,
                                                                       1480, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4612, 0, 3, 4144,
                                                                       1264, 4162, 138, 144,
                                                                       1498, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4648, 0, 3, 4162,
                                                                       1273, 4180, 144, 150,
                                                                       1516, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4684, 0, 3, 4180,
                                                                       1282, 4198, 150, 156,
                                                                       1534, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4720, 0, 3, 4198,
                                                                       1291, 4216, 156, 162,
                                                                       1552, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4756, 0, 3, 4216,
                                                                       1300, 4234, 162, 168,
                                                                       1570, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4792, 0, 3, 4234,
                                                                       1309, 4252, 168, 174,
                                                                       1588, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4828, 0, 3, 4270,
                                                                       1327, 4288, 186, 192,
                                                                       1606, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4864, 0, 3, 4288,
                                                                       1336, 4306, 192, 198,
                                                                       1624, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4900, 0, 3, 4306,
                                                                       1345, 4324, 198, 204,
                                                                       1642, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4936, 0, 3, 4324,
                                                                       1354, 4342, 204, 210,
                                                                       1660, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4972, 0, 3, 4342,
                                                                       1363, 4360, 210, 216,
                                                                       1678, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5008, 0, 3, 4360,
                                                                       1372, 4378, 216, 222,
                                                                       1696, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5044, 0, 3, 4378,
                                                                       1381, 4396, 222, 228,
                                                                       1714, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5080, 0, 3, 4396,
                                                                       1390, 4414, 228, 234,
                                                                       1732, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5116, 0, 3, 4414,
                                                                       1399, 4432, 234, 240,
                                                                       1750, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5152, 0, 3, 4432,
                                                                       1408, 4450, 240, 246,
                                                                       1768, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5188, 0, 3, 4468,
                                                                       1426, 4504, 258, 268,
                                                                       1786, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5248, 0, 3, 4504,
                                                                       1444, 4540, 268, 278,
                                                                       1816, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5308, 0, 3, 4540,
                                                                       1462, 4576, 278, 288,
                                                                       1846, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5368, 0, 3, 4576,
                                                                       1480, 4612, 288, 298,
                                                                       1876, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5428, 0, 3, 4612,
                                                                       1498, 4648, 298, 308,
                                                                       1906, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5488, 0, 3, 4648,
                                                                       1516, 4684, 308, 318,
                                                                       1936, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5548, 0, 3, 4684,
                                                                       1534, 4720, 318, 328,
                                                                       1966, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5608, 0, 3, 4720,
                                                                       1552, 4756, 328, 338,
                                                                       1996, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5668, 0, 3, 4756,
                                                                       1570, 4792, 338, 348,
                                                                       2026, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5728, 0, 3, 4828,
                                                                       1606, 4864, 368, 378,
                                                                       2056, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5788, 0, 3, 4864,
                                                                       1624, 4900, 378, 388,
                                                                       2086, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5848, 0, 3, 4900,
                                                                       1642, 4936, 388, 398,
                                                                       2116, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5908, 0, 3, 4936,
                                                                       1660, 4972, 398, 408,
                                                                       2146, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5968, 0, 3, 4972,
                                                                       1678, 5008, 408, 418,
                                                                       2176, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6028, 0, 3, 5008,
                                                                       1696, 5044, 418, 428,
                                                                       2206, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6088, 0, 3, 5044,
                                                                       1714, 5080, 428, 438,
                                                                       2236, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6148, 0, 3, 5080,
                                                                       1732, 5116, 438, 448,
                                                                       2266, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6208, 0, 3, 5116,
                                                                       1750, 5152, 448, 458,
                                                                       2296, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6268, 0, 3, 5188,
                                                                       1786, 5248, 478, 493,
                                                                       2326, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6358, 0, 3, 5248,
                                                                       1816, 5308, 493, 508,
                                                                       2371, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6448, 0, 3, 5308,
                                                                       1846, 5368, 508, 523,
                                                                       2416, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6538, 0, 3, 5368,
                                                                       1876, 5428, 523, 538,
                                                                       2461, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6628, 0, 3, 5428,
                                                                       1906, 5488, 538, 553,
                                                                       2506, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6718, 0, 3, 5488,
                                                                       1936, 5548, 553, 568,
                                                                       2551, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6808, 0, 3, 5548,
                                                                       1966, 5608, 568, 583,
                                                                       2596, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6898, 0, 3, 5608,
                                                                       1996, 5668, 583, 598,
                                                                       2641, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6988, 0, 3, 5728,
                                                                       2056, 5788, 628, 643,
                                                                       2686, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7078, 0, 3, 5788,
                                                                       2086, 5848, 643, 658,
                                                                       2731, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7168, 0, 3, 5848,
                                                                       2116, 5908, 658, 673,
                                                                       2776, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7258, 0, 3, 5908,
                                                                       2146, 5968, 673, 688,
                                                                       2821, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7348, 0, 3, 5968,
                                                                       2176, 6028, 688, 703,
                                                                       2866, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7438, 0, 3, 6028,
                                                                       2206, 6088, 703, 718,
                                                                       2911, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7528, 0, 3, 6088,
                                                                       2236, 6148, 718, 733,
                                                                       2956, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7618, 0, 3, 6148,
                                                                       2266, 6208, 733, 748,
                                                                       3001, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7708, 0, 3, 6268,
                                                                       2326, 6358, 778, 799,
                                                                       3046, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7834, 0, 3, 6358,
                                                                       2371, 6448, 799, 820,
                                                                       3109, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7960, 0, 3, 6448,
                                                                       2416, 6538, 820, 841,
                                                                       3172, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8086, 0, 3, 6538,
                                                                       2461, 6628, 841, 862,
                                                                       3235, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8212, 0, 3, 6628,
                                                                       2506, 6718, 862, 883,
                                                                       3298, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8338, 0, 3, 6718,
                                                                       2551, 6808, 883, 904,
                                                                       3361, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8464, 0, 3, 6808,
                                                                       2596, 6898, 904, 925,
                                                                       3424, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8590, 0, 3, 6988,
                                                                       2686, 7078, 967, 988,
                                                                       3487, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8716, 0, 3, 7078,
                                                                       2731, 7168, 988, 1009,
                                                                       3550, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8842, 0, 3, 7168,
                                                                       2776, 7258, 1009, 1030,
                                                                       3613, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8968, 0, 3, 7258,
                                                                       2821, 7348, 1030, 1051,
                                                                       3676, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9094, 0, 3, 7348,
                                                                       2866, 7438, 1051, 1072,
                                                                       3739, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9220, 0, 3, 7438,
                                                                       2911, 7528, 1072, 1093,
                                                                       3802, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9346, 0, 3, 7528,
                                                                       2956, 7618, 1093, 1114,
                                                                       3865, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9472, 3, 1156,
                                                                       1159, 3940, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9482, 3, 1159,
                                                                       1162, 3946, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9492, 3, 1162,
                                                                       1165, 3952, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9502, 3, 1165,
                                                                       1168, 3958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9512, 3, 1168,
                                                                       1171, 3964, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9522, 3, 1171,
                                                                       1174, 3970, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9532, 3, 1174,
                                                                       1177, 3976, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9542, 3, 1177,
                                                                       1180, 3982, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9552, 3, 1180,
                                                                       1183, 3988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9562, 3, 1183,
                                                                       1186, 3994, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9572, 3, 1192,
                                                                       1195, 4012, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9582, 3, 1195,
                                                                       1198, 4018, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9592, 3, 1198,
                                                                       1201, 4024, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9602, 3, 1201,
                                                                       1204, 4030, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9612, 3, 1204,
                                                                       1207, 4036, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9622, 3, 1207,
                                                                       1210, 4042, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9632, 3, 1210,
                                                                       1213, 4048, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9642, 3, 1213,
                                                                       1216, 4054, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9652, 3, 1216,
                                                                       1219, 4060, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9662, 3, 1219,
                                                                       1222, 4066, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9672, 0, 3, 9472,
                                                                       3940, 9482, 4108, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9702, 0, 3, 9482,
                                                                       3946, 9492, 4126, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9732, 0, 3, 9492,
                                                                       3952, 9502, 4144, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9762, 0, 3, 9502,
                                                                       3958, 9512, 4162, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9792, 0, 3, 9512,
                                                                       3964, 9522, 4180, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9822, 0, 3, 9522,
                                                                       3970, 9532, 4198, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9852, 0, 3, 9532,
                                                                       3976, 9542, 4216, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9882, 0, 3, 9542,
                                                                       3982, 9552, 4234, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9912, 0, 3, 9552,
                                                                       3988, 9562, 4252, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9942, 0, 3, 9572,
                                                                       4012, 9582, 4306, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9972, 0, 3, 9582,
                                                                       4018, 9592, 4324, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10002, 0, 3, 9592,
                                                                       4024, 9602, 4342, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10032, 0, 3, 9602,
                                                                       4030, 9612, 4360, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10062, 0, 3, 9612,
                                                                       4036, 9622, 4378, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10092, 0, 3, 9622,
                                                                       4042, 9632, 4396, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10122, 0, 3, 9632,
                                                                       4048, 9642, 4414, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10152, 0, 3, 9642,
                                                                       4054, 9652, 4432, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10182, 0, 3, 9652,
                                                                       4060, 9662, 4450, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10212, 0, 3, 9672,
                                                                       4108, 9702, 1426, 1444,
                                                                       4540, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10272, 0, 3, 9702,
                                                                       4126, 9732, 1444, 1462,
                                                                       4576, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10332, 0, 3, 9732,
                                                                       4144, 9762, 1462, 1480,
                                                                       4612, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10392, 0, 3, 9762,
                                                                       4162, 9792, 1480, 1498,
                                                                       4648, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10452, 0, 3, 9792,
                                                                       4180, 9822, 1498, 1516,
                                                                       4684, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10512, 0, 3, 9822,
                                                                       4198, 9852, 1516, 1534,
                                                                       4720, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10572, 0, 3, 9852,
                                                                       4216, 9882, 1534, 1552,
                                                                       4756, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10632, 0, 3, 9882,
                                                                       4234, 9912, 1552, 1570,
                                                                       4792, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10692, 0, 3, 9942,
                                                                       4306, 9972, 1606, 1624,
                                                                       4900, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10752, 0, 3, 9972,
                                                                       4324, 10002, 1624, 1642,
                                                                       4936, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10812, 0, 3,
                                                                       10002, 4342, 10032, 1642,
                                                                       1660, 4972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10872, 0, 3,
                                                                       10032, 4360, 10062, 1660,
                                                                       1678, 5008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10932, 0, 3,
                                                                       10062, 4378, 10092, 1678,
                                                                       1696, 5044, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10992, 0, 3,
                                                                       10092, 4396, 10122, 1696,
                                                                       1714, 5080, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11052, 0, 3,
                                                                       10122, 4414, 10152, 1714,
                                                                       1732, 5116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11112, 0, 3,
                                                                       10152, 4432, 10182, 1732,
                                                                       1750, 5152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11172, 0, 3,
                                                                       10212, 4540, 10272, 1786,
                                                                       1816, 5308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11272, 0, 3,
                                                                       10272, 4576, 10332, 1816,
                                                                       1846, 5368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11372, 0, 3,
                                                                       10332, 4612, 10392, 1846,
                                                                       1876, 5428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11472, 0, 3,
                                                                       10392, 4648, 10452, 1876,
                                                                       1906, 5488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11572, 0, 3,
                                                                       10452, 4684, 10512, 1906,
                                                                       1936, 5548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11672, 0, 3,
                                                                       10512, 4720, 10572, 1936,
                                                                       1966, 5608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11772, 0, 3,
                                                                       10572, 4756, 10632, 1966,
                                                                       1996, 5668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11872, 0, 3,
                                                                       10692, 4900, 10752, 2056,
                                                                       2086, 5848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11972, 0, 3,
                                                                       10752, 4936, 10812, 2086,
                                                                       2116, 5908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12072, 0, 3,
                                                                       10812, 4972, 10872, 2116,
                                                                       2146, 5968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12172, 0, 3,
                                                                       10872, 5008, 10932, 2146,
                                                                       2176, 6028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12272, 0, 3,
                                                                       10932, 5044, 10992, 2176,
                                                                       2206, 6088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12372, 0, 3,
                                                                       10992, 5080, 11052, 2206,
                                                                       2236, 6148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12472, 0, 3,
                                                                       11052, 5116, 11112, 2236,
                                                                       2266, 6208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12572, 0, 3,
                                                                       11172, 5308, 11272, 2326,
                                                                       2371, 6448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12722, 0, 3,
                                                                       11272, 5368, 11372, 2371,
                                                                       2416, 6538, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12872, 0, 3,
                                                                       11372, 5428, 11472, 2416,
                                                                       2461, 6628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13022, 0, 3,
                                                                       11472, 5488, 11572, 2461,
                                                                       2506, 6718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13172, 0, 3,
                                                                       11572, 5548, 11672, 2506,
                                                                       2551, 6808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13322, 0, 3,
                                                                       11672, 5608, 11772, 2551,
                                                                       2596, 6898, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13472, 0, 3,
                                                                       11872, 5848, 11972, 2686,
                                                                       2731, 7168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13622, 0, 3,
                                                                       11972, 5908, 12072, 2731,
                                                                       2776, 7258, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13772, 0, 3,
                                                                       12072, 5968, 12172, 2776,
                                                                       2821, 7348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13922, 0, 3,
                                                                       12172, 6028, 12272, 2821,
                                                                       2866, 7438, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14072, 0, 3,
                                                                       12272, 6088, 12372, 2866,
                                                                       2911, 7528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14222, 0, 3,
                                                                       12372, 6148, 12472, 2911,
                                                                       2956, 7618, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14372, 0, 3,
                                                                       12572, 6448, 12722, 3046,
                                                                       3109, 7960, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14582, 0, 3,
                                                                       12722, 6538, 12872, 3109,
                                                                       3172, 8086, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14792, 0, 3,
                                                                       12872, 6628, 13022, 3172,
                                                                       3235, 8212, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15002, 0, 3,
                                                                       13022, 6718, 13172, 3235,
                                                                       3298, 8338, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15212, 0, 3,
                                                                       13172, 6808, 13322, 3298,
                                                                       3361, 8464, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15422, 0, 3,
                                                                       13472, 7168, 13622, 3487,
                                                                       3550, 8842, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15632, 0, 3,
                                                                       13622, 7258, 13772, 3550,
                                                                       3613, 8968, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15842, 0, 3,
                                                                       13772, 7348, 13922, 3613,
                                                                       3676, 9094, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16052, 0, 3,
                                                                       13922, 7438, 14072, 3676,
                                                                       3739, 9220, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16262, 0, 3,
                                                                       14072, 7528, 14222, 3739,
                                                                       3802, 9346, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16472, 3, 3928,
                                                                       3934, 9472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16487, 3, 3934,
                                                                       3940, 9482, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16502, 3, 3940,
                                                                       3946, 9492, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16517, 3, 3946,
                                                                       3952, 9502, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16532, 3, 3952,
                                                                       3958, 9512, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16547, 3, 3958,
                                                                       3964, 9522, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16562, 3, 3964,
                                                                       3970, 9532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16577, 3, 3970,
                                                                       3976, 9542, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16592, 3, 3976,
                                                                       3982, 9552, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16607, 3, 3982,
                                                                       3988, 9562, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16622, 3, 4000,
                                                                       4006, 9572, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16637, 3, 4006,
                                                                       4012, 9582, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16652, 3, 4012,
                                                                       4018, 9592, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16667, 3, 4018,
                                                                       4024, 9602, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16682, 3, 4024,
                                                                       4030, 9612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16697, 3, 4030,
                                                                       4036, 9622, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16712, 3, 4036,
                                                                       4042, 9632, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16727, 3, 4042,
                                                                       4048, 9642, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16742, 3, 4048,
                                                                       4054, 9652, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16757, 3, 4054,
                                                                       4060, 9662, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16772, 0, 3,
                                                                       16472, 9472, 16487, 4072,
                                                                       4090, 9672, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16817, 0, 3,
                                                                       16487, 9482, 16502, 4090,
                                                                       4108, 9702, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16862, 0, 3,
                                                                       16502, 9492, 16517, 4108,
                                                                       4126, 9732, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16907, 0, 3,
                                                                       16517, 9502, 16532, 4126,
                                                                       4144, 9762, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16952, 0, 3,
                                                                       16532, 9512, 16547, 4144,
                                                                       4162, 9792, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16997, 0, 3,
                                                                       16547, 9522, 16562, 4162,
                                                                       4180, 9822, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17042, 0, 3,
                                                                       16562, 9532, 16577, 4180,
                                                                       4198, 9852, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17087, 0, 3,
                                                                       16577, 9542, 16592, 4198,
                                                                       4216, 9882, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17132, 0, 3,
                                                                       16592, 9552, 16607, 4216,
                                                                       4234, 9912, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17177, 0, 3,
                                                                       16622, 9572, 16637, 4270,
                                                                       4288, 9942, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17222, 0, 3,
                                                                       16637, 9582, 16652, 4288,
                                                                       4306, 9972, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17267, 0, 3,
                                                                       16652, 9592, 16667, 4306,
                                                                       4324, 10002, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17312, 0, 3,
                                                                       16667, 9602, 16682, 4324,
                                                                       4342, 10032, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17357, 0, 3,
                                                                       16682, 9612, 16697, 4342,
                                                                       4360, 10062, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17402, 0, 3,
                                                                       16697, 9622, 16712, 4360,
                                                                       4378, 10092, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17447, 0, 3,
                                                                       16712, 9632, 16727, 4378,
                                                                       4396, 10122, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17492, 0, 3,
                                                                       16727, 9642, 16742, 4396,
                                                                       4414, 10152, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17537, 0, 3,
                                                                       16742, 9652, 16757, 4414,
                                                                       4432, 10182, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17582, 0, 3,
                                                                       16772, 9672, 16817, 4468,
                                                                       4504, 10212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17672, 0, 3,
                                                                       16817, 9702, 16862, 4504,
                                                                       4540, 10272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17762, 0, 3,
                                                                       16862, 9732, 16907, 4540,
                                                                       4576, 10332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17852, 0, 3,
                                                                       16907, 9762, 16952, 4576,
                                                                       4612, 10392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17942, 0, 3,
                                                                       16952, 9792, 16997, 4612,
                                                                       4648, 10452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18032, 0, 3,
                                                                       16997, 9822, 17042, 4648,
                                                                       4684, 10512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18122, 0, 3,
                                                                       17042, 9852, 17087, 4684,
                                                                       4720, 10572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18212, 0, 3,
                                                                       17087, 9882, 17132, 4720,
                                                                       4756, 10632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18302, 0, 3,
                                                                       17177, 9942, 17222, 4828,
                                                                       4864, 10692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18392, 0, 3,
                                                                       17222, 9972, 17267, 4864,
                                                                       4900, 10752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18482, 0, 3,
                                                                       17267, 10002, 17312, 4900,
                                                                       4936, 10812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18572, 0, 3,
                                                                       17312, 10032, 17357, 4936,
                                                                       4972, 10872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18662, 0, 3,
                                                                       17357, 10062, 17402, 4972,
                                                                       5008, 10932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18752, 0, 3,
                                                                       17402, 10092, 17447, 5008,
                                                                       5044, 10992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18842, 0, 3,
                                                                       17447, 10122, 17492, 5044,
                                                                       5080, 11052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18932, 0, 3,
                                                                       17492, 10152, 17537, 5080,
                                                                       5116, 11112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19022, 0, 3,
                                                                       17582, 10212, 17672, 5188,
                                                                       5248, 11172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19172, 0, 3,
                                                                       17672, 10272, 17762, 5248,
                                                                       5308, 11272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19322, 0, 3,
                                                                       17762, 10332, 17852, 5308,
                                                                       5368, 11372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19472, 0, 3,
                                                                       17852, 10392, 17942, 5368,
                                                                       5428, 11472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19622, 0, 3,
                                                                       17942, 10452, 18032, 5428,
                                                                       5488, 11572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19772, 0, 3,
                                                                       18032, 10512, 18122, 5488,
                                                                       5548, 11672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19922, 0, 3,
                                                                       18122, 10572, 18212, 5548,
                                                                       5608, 11772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20072, 0, 3,
                                                                       18302, 10692, 18392, 5728,
                                                                       5788, 11872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20222, 0, 3,
                                                                       18392, 10752, 18482, 5788,
                                                                       5848, 11972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20372, 0, 3,
                                                                       18482, 10812, 18572, 5848,
                                                                       5908, 12072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20522, 0, 3,
                                                                       18572, 10872, 18662, 5908,
                                                                       5968, 12172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20672, 0, 3,
                                                                       18662, 10932, 18752, 5968,
                                                                       6028, 12272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20822, 0, 3,
                                                                       18752, 10992, 18842, 6028,
                                                                       6088, 12372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20972, 0, 3,
                                                                       18842, 11052, 18932, 6088,
                                                                       6148, 12472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21122, 0, 3,
                                                                       19022, 11172, 19172, 6268,
                                                                       6358, 12572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21347, 0, 3,
                                                                       19172, 11272, 19322, 6358,
                                                                       6448, 12722, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21572, 0, 3,
                                                                       19322, 11372, 19472, 6448,
                                                                       6538, 12872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21797, 0, 3,
                                                                       19472, 11472, 19622, 6538,
                                                                       6628, 13022, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22022, 0, 3,
                                                                       19622, 11572, 19772, 6628,
                                                                       6718, 13172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22247, 0, 3,
                                                                       19772, 11672, 19922, 6718,
                                                                       6808, 13322, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22472, 0, 3,
                                                                       20072, 11872, 20222, 6988,
                                                                       7078, 13472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22697, 0, 3,
                                                                       20222, 11972, 20372, 7078,
                                                                       7168, 13622, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22922, 0, 3,
                                                                       20372, 12072, 20522, 7168,
                                                                       7258, 13772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23147, 0, 3,
                                                                       20522, 12172, 20672, 7258,
                                                                       7348, 13922, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23372, 0, 3,
                                                                       20672, 12272, 20822, 7348,
                                                                       7438, 14072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23597, 0, 3,
                                                                       20822, 12372, 20972, 7438,
                                                                       7528, 14222, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 23822, 0, 3,
                                                                       21122, 12572, 21347, 7708,
                                                                       7834, 14372, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24137, 0, 3,
                                                                       21347, 12722, 21572, 7834,
                                                                       7960, 14582, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24452, 0, 3,
                                                                       21572, 12872, 21797, 7960,
                                                                       8086, 14792, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24767, 0, 3,
                                                                       21797, 13022, 22022, 8086,
                                                                       8212, 15002, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25082, 0, 3,
                                                                       22022, 13172, 22247, 8212,
                                                                       8338, 15212, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25397, 0, 3,
                                                                       22472, 13472, 22697, 8590,
                                                                       8716, 15422, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25712, 0, 3,
                                                                       22697, 13622, 22922, 8716,
                                                                       8842, 15632, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 26027, 0, 3,
                                                                       22922, 13772, 23147, 8842,
                                                                       8968, 15842, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 26342, 0, 3,
                                                                       23147, 13922, 23372, 8968,
                                                                       9094, 16052, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 26657, 0, 3,
                                                                       23372, 14072, 23597, 9094,
                                                                       9220, 16262, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26972, 3, 9472,
                                                                       9482, 16502, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26993, 3, 9482,
                                                                       9492, 16517, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27014, 3, 9492,
                                                                       9502, 16532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27035, 3, 9502,
                                                                       9512, 16547, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27056, 3, 9512,
                                                                       9522, 16562, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27077, 3, 9522,
                                                                       9532, 16577, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27098, 3, 9532,
                                                                       9542, 16592, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27119, 3, 9542,
                                                                       9552, 16607, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27140, 3, 9572,
                                                                       9582, 16652, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27161, 3, 9582,
                                                                       9592, 16667, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27182, 3, 9592,
                                                                       9602, 16682, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27203, 3, 9602,
                                                                       9612, 16697, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27224, 3, 9612,
                                                                       9622, 16712, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27245, 3, 9622,
                                                                       9632, 16727, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27266, 3, 9632,
                                                                       9642, 16742, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 27287, 3, 9642,
                                                                       9652, 16757, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27308, 0, 3,
                                                                       26972, 16502, 26993, 9672,
                                                                       9702, 16862, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27371, 0, 3,
                                                                       26993, 16517, 27014, 9702,
                                                                       9732, 16907, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27434, 0, 3,
                                                                       27014, 16532, 27035, 9732,
                                                                       9762, 16952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27497, 0, 3,
                                                                       27035, 16547, 27056, 9762,
                                                                       9792, 16997, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27560, 0, 3,
                                                                       27056, 16562, 27077, 9792,
                                                                       9822, 17042, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27623, 0, 3,
                                                                       27077, 16577, 27098, 9822,
                                                                       9852, 17087, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27686, 0, 3,
                                                                       27098, 16592, 27119, 9852,
                                                                       9882, 17132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27749, 0, 3,
                                                                       27140, 16652, 27161, 9942,
                                                                       9972, 17267, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27812, 0, 3,
                                                                       27161, 16667, 27182, 9972,
                                                                       10002, 17312, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27875, 0, 3,
                                                                       27182, 16682, 27203,
                                                                       10002, 10032, 17357,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 27938, 0, 3,
                                                                       27203, 16697, 27224,
                                                                       10032, 10062, 17402,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 28001, 0, 3,
                                                                       27224, 16712, 27245,
                                                                       10062, 10092, 17447,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 28064, 0, 3,
                                                                       27245, 16727, 27266,
                                                                       10092, 10122, 17492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 28127, 0, 3,
                                                                       27266, 16742, 27287,
                                                                       10122, 10152, 17537,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 28190, 0, 3,
                                                                       27308, 16862, 27371,
                                                                       10212, 10272, 17762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 28316, 0, 3,
                                                                       27371, 16907, 27434,
                                                                       10272, 10332, 17852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 28442, 0, 3,
                                                                       27434, 16952, 27497,
                                                                       10332, 10392, 17942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 28568, 0, 3,
                                                                       27497, 16997, 27560,
                                                                       10392, 10452, 18032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 28694, 0, 3,
                                                                       27560, 17042, 27623,
                                                                       10452, 10512, 18122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 28820, 0, 3,
                                                                       27623, 17087, 27686,
                                                                       10512, 10572, 18212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 28946, 0, 3,
                                                                       27749, 17267, 27812,
                                                                       10692, 10752, 18482,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 29072, 0, 3,
                                                                       27812, 17312, 27875,
                                                                       10752, 10812, 18572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 29198, 0, 3,
                                                                       27875, 17357, 27938,
                                                                       10812, 10872, 18662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 29324, 0, 3,
                                                                       27938, 17402, 28001,
                                                                       10872, 10932, 18752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 29450, 0, 3,
                                                                       28001, 17447, 28064,
                                                                       10932, 10992, 18842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 29576, 0, 3,
                                                                       28064, 17492, 28127,
                                                                       10992, 11052, 18932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 29702, 0, 3,
                                                                       28190, 17762, 28316,
                                                                       11172, 11272, 19322,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 29912, 0, 3,
                                                                       28316, 17852, 28442,
                                                                       11272, 11372, 19472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 30122, 0, 3,
                                                                       28442, 17942, 28568,
                                                                       11372, 11472, 19622,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 30332, 0, 3,
                                                                       28568, 18032, 28694,
                                                                       11472, 11572, 19772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 30542, 0, 3,
                                                                       28694, 18122, 28820,
                                                                       11572, 11672, 19922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 30752, 0, 3,
                                                                       28946, 18482, 29072,
                                                                       11872, 11972, 20372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 30962, 0, 3,
                                                                       29072, 18572, 29198,
                                                                       11972, 12072, 20522,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 31172, 0, 3,
                                                                       29198, 18662, 29324,
                                                                       12072, 12172, 20672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 31382, 0, 3,
                                                                       29324, 18752, 29450,
                                                                       12172, 12272, 20822,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 31592, 0, 3,
                                                                       29450, 18842, 29576,
                                                                       12272, 12372, 20972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 31802, 0, 3,
                                                                       29702, 19322, 29912,
                                                                       12572, 12722, 21572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 32117, 0, 3,
                                                                       29912, 19472, 30122,
                                                                       12722, 12872, 21797,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 32432, 0, 3,
                                                                       30122, 19622, 30332,
                                                                       12872, 13022, 22022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 32747, 0, 3,
                                                                       30332, 19772, 30542,
                                                                       13022, 13172, 22247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 33062, 0, 3,
                                                                       30752, 20372, 30962,
                                                                       13472, 13622, 22922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 33377, 0, 3,
                                                                       30962, 20522, 31172,
                                                                       13622, 13772, 23147,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 33692, 0, 3,
                                                                       31172, 20672, 31382,
                                                                       13772, 13922, 23372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 34007, 0, 3,
                                                                       31382, 20822, 31592,
                                                                       13922, 14072, 23597,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 34322, 0, 3,
                                                                       31802, 21572, 32117,
                                                                       14372, 14582, 24452,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 34763, 0, 3,
                                                                       32117, 21797, 32432,
                                                                       14582, 14792, 24767,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 35204, 0, 3,
                                                                       32432, 22022, 32747,
                                                                       14792, 15002, 25082,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 35645, 0, 3,
                                                                       33062, 22922, 33377,
                                                                       15422, 15632, 26027,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 36086, 0, 3,
                                                                       33377, 23147, 33692,
                                                                       15632, 15842, 26342,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 36527, 0, 3,
                                                                       33692, 23372, 34007,
                                                                       15842, 16052, 26657,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 36968, 3, 16472,
                                                                       16487, 26972, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 36996, 3, 16487,
                                                                       16502, 26993, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37024, 3, 16502,
                                                                       16517, 27014, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37052, 3, 16517,
                                                                       16532, 27035, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37080, 3, 16532,
                                                                       16547, 27056, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37108, 3, 16547,
                                                                       16562, 27077, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37136, 3, 16562,
                                                                       16577, 27098, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37164, 3, 16577,
                                                                       16592, 27119, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37192, 3, 16622,
                                                                       16637, 27140, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37220, 3, 16637,
                                                                       16652, 27161, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37248, 3, 16652,
                                                                       16667, 27182, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37276, 3, 16667,
                                                                       16682, 27203, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37304, 3, 16682,
                                                                       16697, 27224, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37332, 3, 16697,
                                                                       16712, 27245, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37360, 3, 16712,
                                                                       16727, 27266, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37388, 3, 16727,
                                                                       16742, 27287, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 37416, 0, 3,
                                                                       36968, 26972, 36996,
                                                                       16772, 16817, 27308,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 37500, 0, 3,
                                                                       36996, 26993, 37024,
                                                                       16817, 16862, 27371,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 37584, 0, 3,
                                                                       37024, 27014, 37052,
                                                                       16862, 16907, 27434,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 37668, 0, 3,
                                                                       37052, 27035, 37080,
                                                                       16907, 16952, 27497,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 37752, 0, 3,
                                                                       37080, 27056, 37108,
                                                                       16952, 16997, 27560,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 37836, 0, 3,
                                                                       37108, 27077, 37136,
                                                                       16997, 17042, 27623,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 37920, 0, 3,
                                                                       37136, 27098, 37164,
                                                                       17042, 17087, 27686,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38004, 0, 3,
                                                                       37192, 27140, 37220,
                                                                       17177, 17222, 27749,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38088, 0, 3,
                                                                       37220, 27161, 37248,
                                                                       17222, 17267, 27812,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38172, 0, 3,
                                                                       37248, 27182, 37276,
                                                                       17267, 17312, 27875,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38256, 0, 3,
                                                                       37276, 27203, 37304,
                                                                       17312, 17357, 27938,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38340, 0, 3,
                                                                       37304, 27224, 37332,
                                                                       17357, 17402, 28001,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38424, 0, 3,
                                                                       37332, 27245, 37360,
                                                                       17402, 17447, 28064,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38508, 0, 3,
                                                                       37360, 27266, 37388,
                                                                       17447, 17492, 28127,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 38592, 0, 3,
                                                                       37416, 27308, 37500,
                                                                       17582, 17672, 28190,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 38760, 0, 3,
                                                                       37500, 27371, 37584,
                                                                       17672, 17762, 28316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 38928, 0, 3,
                                                                       37584, 27434, 37668,
                                                                       17762, 17852, 28442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 39096, 0, 3,
                                                                       37668, 27497, 37752,
                                                                       17852, 17942, 28568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 39264, 0, 3,
                                                                       37752, 27560, 37836,
                                                                       17942, 18032, 28694,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 39432, 0, 3,
                                                                       37836, 27623, 37920,
                                                                       18032, 18122, 28820,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 39600, 0, 3,
                                                                       38004, 27749, 38088,
                                                                       18302, 18392, 28946,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 39768, 0, 3,
                                                                       38088, 27812, 38172,
                                                                       18392, 18482, 29072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 39936, 0, 3,
                                                                       38172, 27875, 38256,
                                                                       18482, 18572, 29198,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 40104, 0, 3,
                                                                       38256, 27938, 38340,
                                                                       18572, 18662, 29324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 40272, 0, 3,
                                                                       38340, 28001, 38424,
                                                                       18662, 18752, 29450,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 40440, 0, 3,
                                                                       38424, 28064, 38508,
                                                                       18752, 18842, 29576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 40608, 0, 3,
                                                                       38592, 28190, 38760,
                                                                       19022, 19172, 29702,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 40888, 0, 3,
                                                                       38760, 28316, 38928,
                                                                       19172, 19322, 29912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 41168, 0, 3,
                                                                       38928, 28442, 39096,
                                                                       19322, 19472, 30122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 41448, 0, 3,
                                                                       39096, 28568, 39264,
                                                                       19472, 19622, 30332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 41728, 0, 3,
                                                                       39264, 28694, 39432,
                                                                       19622, 19772, 30542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 42008, 0, 3,
                                                                       39600, 28946, 39768,
                                                                       20072, 20222, 30752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 42288, 0, 3,
                                                                       39768, 29072, 39936,
                                                                       20222, 20372, 30962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 42568, 0, 3,
                                                                       39936, 29198, 40104,
                                                                       20372, 20522, 31172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 42848, 0, 3,
                                                                       40104, 29324, 40272,
                                                                       20522, 20672, 31382,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 43128, 0, 3,
                                                                       40272, 29450, 40440,
                                                                       20672, 20822, 31592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 43408, 0, 3,
                                                                       40608, 29702, 40888,
                                                                       21122, 21347, 31802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 43828, 0, 3,
                                                                       40888, 29912, 41168,
                                                                       21347, 21572, 32117,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 44248, 0, 3,
                                                                       41168, 30122, 41448,
                                                                       21572, 21797, 32432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 44668, 0, 3,
                                                                       41448, 30332, 41728,
                                                                       21797, 22022, 32747,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 45088, 0, 3,
                                                                       42008, 30752, 42288,
                                                                       22472, 22697, 33062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 45508, 0, 3,
                                                                       42288, 30962, 42568,
                                                                       22697, 22922, 33377,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 45928, 0, 3,
                                                                       42568, 31172, 42848,
                                                                       22922, 23147, 33692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 46348, 0, 3,
                                                                       42848, 31382, 43128,
                                                                       23147, 23372, 34007,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 46768, 0, 3,
                                                                       43408, 31802, 43828,
                                                                       23822, 24137, 34322,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 47356, 0, 3,
                                                                       43828, 32117, 44248,
                                                                       24137, 24452, 34763,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 47944, 0, 3,
                                                                       44248, 32432, 44668,
                                                                       24452, 24767, 35204,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 48532, 0, 3,
                                                                       45088, 33062, 45508,
                                                                       25397, 25712, 35645,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 49120, 0, 3,
                                                                       45508, 33377, 45928,
                                                                       25712, 26027, 36086,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 49708, 0, 3,
                                                                       45928, 33692, 46348,
                                                                       26027, 26342, 36527,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50296, 3, 26972,
                                                                       26993, 37024, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50332, 3, 26993,
                                                                       27014, 37052, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50368, 3, 27014,
                                                                       27035, 37080, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50404, 3, 27035,
                                                                       27056, 37108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50440, 3, 27056,
                                                                       27077, 37136, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50476, 3, 27077,
                                                                       27098, 37164, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50512, 3, 27140,
                                                                       27161, 37248, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50548, 3, 27161,
                                                                       27182, 37276, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50584, 3, 27182,
                                                                       27203, 37304, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50620, 3, 27203,
                                                                       27224, 37332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50656, 3, 27224,
                                                                       27245, 37360, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 50692, 3, 27245,
                                                                       27266, 37388, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 50728, 0, 3,
                                                                       50296, 37024, 50332,
                                                                       27308, 27371, 37584,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 50836, 0, 3,
                                                                       50332, 37052, 50368,
                                                                       27371, 27434, 37668,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 50944, 0, 3,
                                                                       50368, 37080, 50404,
                                                                       27434, 27497, 37752,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 51052, 0, 3,
                                                                       50404, 37108, 50440,
                                                                       27497, 27560, 37836,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 51160, 0, 3,
                                                                       50440, 37136, 50476,
                                                                       27560, 27623, 37920,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 51268, 0, 3,
                                                                       50512, 37248, 50548,
                                                                       27749, 27812, 38172,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 51376, 0, 3,
                                                                       50548, 37276, 50584,
                                                                       27812, 27875, 38256,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 51484, 0, 3,
                                                                       50584, 37304, 50620,
                                                                       27875, 27938, 38340,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 51592, 0, 3,
                                                                       50620, 37332, 50656,
                                                                       27938, 28001, 38424,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 51700, 0, 3,
                                                                       50656, 37360, 50692,
                                                                       28001, 28064, 38508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 51808, 0, 3,
                                                                       50728, 37584, 50836,
                                                                       28190, 28316, 38928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 52024, 0, 3,
                                                                       50836, 37668, 50944,
                                                                       28316, 28442, 39096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 52240, 0, 3,
                                                                       50944, 37752, 51052,
                                                                       28442, 28568, 39264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 52456, 0, 3,
                                                                       51052, 37836, 51160,
                                                                       28568, 28694, 39432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 52672, 0, 3,
                                                                       51268, 38172, 51376,
                                                                       28946, 29072, 39936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 52888, 0, 3,
                                                                       51376, 38256, 51484,
                                                                       29072, 29198, 40104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 53104, 0, 3,
                                                                       51484, 38340, 51592,
                                                                       29198, 29324, 40272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 53320, 0, 3,
                                                                       51592, 38424, 51700,
                                                                       29324, 29450, 40440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 53536, 0, 3,
                                                                       51808, 38928, 52024,
                                                                       29702, 29912, 41168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 53896, 0, 3,
                                                                       52024, 39096, 52240,
                                                                       29912, 30122, 41448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 54256, 0, 3,
                                                                       52240, 39264, 52456,
                                                                       30122, 30332, 41728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 54616, 0, 3,
                                                                       52672, 39936, 52888,
                                                                       30752, 30962, 42568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 54976, 0, 3,
                                                                       52888, 40104, 53104,
                                                                       30962, 31172, 42848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 55336, 0, 3,
                                                                       53104, 40272, 53320,
                                                                       31172, 31382, 43128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 55696, 0, 3,
                                                                       53536, 41168, 53896,
                                                                       31802, 32117, 44248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 56236, 0, 3,
                                                                       53896, 41448, 54256,
                                                                       32117, 32432, 44668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 56776, 0, 3,
                                                                       54616, 42568, 54976,
                                                                       33062, 33377, 45928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 57316, 0, 3,
                                                                       54976, 42848, 55336,
                                                                       33377, 33692, 46348,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 57856, 0, 3,
                                                                       55696, 44248, 56236,
                                                                       34322, 34763, 47944,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 58612, 0, 3,
                                                                       56776, 45928, 57316,
                                                                       35645, 36086, 49708,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59368, 3, 36968,
                                                                       36996, 50296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59413, 3, 36996,
                                                                       37024, 50332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59458, 3, 37024,
                                                                       37052, 50368, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59503, 3, 37052,
                                                                       37080, 50404, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59548, 3, 37080,
                                                                       37108, 50440, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59593, 3, 37108,
                                                                       37136, 50476, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59638, 3, 37192,
                                                                       37220, 50512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59683, 3, 37220,
                                                                       37248, 50548, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59728, 3, 37248,
                                                                       37276, 50584, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59773, 3, 37276,
                                                                       37304, 50620, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59818, 3, 37304,
                                                                       37332, 50656, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 59863, 3, 37332,
                                                                       37360, 50692, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 59908, 0, 3,
                                                                       59368, 50296, 59413,
                                                                       37416, 37500, 50728,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 60043, 0, 3,
                                                                       59413, 50332, 59458,
                                                                       37500, 37584, 50836,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 60178, 0, 3,
                                                                       59458, 50368, 59503,
                                                                       37584, 37668, 50944,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 60313, 0, 3,
                                                                       59503, 50404, 59548,
                                                                       37668, 37752, 51052,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 60448, 0, 3,
                                                                       59548, 50440, 59593,
                                                                       37752, 37836, 51160,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 60583, 0, 3,
                                                                       59638, 50512, 59683,
                                                                       38004, 38088, 51268,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 60718, 0, 3,
                                                                       59683, 50548, 59728,
                                                                       38088, 38172, 51376,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 60853, 0, 3,
                                                                       59728, 50584, 59773,
                                                                       38172, 38256, 51484,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 60988, 0, 3,
                                                                       59773, 50620, 59818,
                                                                       38256, 38340, 51592,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 61123, 0, 3,
                                                                       59818, 50656, 59863,
                                                                       38340, 38424, 51700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 61258, 0, 3,
                                                                       59908, 50728, 60043,
                                                                       38592, 38760, 51808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 61528, 0, 3,
                                                                       60043, 50836, 60178,
                                                                       38760, 38928, 52024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 61798, 0, 3,
                                                                       60178, 50944, 60313,
                                                                       38928, 39096, 52240,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 62068, 0, 3,
                                                                       60313, 51052, 60448,
                                                                       39096, 39264, 52456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 62338, 0, 3,
                                                                       60583, 51268, 60718,
                                                                       39600, 39768, 52672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 62608, 0, 3,
                                                                       60718, 51376, 60853,
                                                                       39768, 39936, 52888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 62878, 0, 3,
                                                                       60853, 51484, 60988,
                                                                       39936, 40104, 53104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 63148, 0, 3,
                                                                       60988, 51592, 61123,
                                                                       40104, 40272, 53320,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 63418, 0, 3,
                                                                       61258, 51808, 61528,
                                                                       40608, 40888, 53536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 63868, 0, 3,
                                                                       61528, 52024, 61798,
                                                                       40888, 41168, 53896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 64318, 0, 3,
                                                                       61798, 52240, 62068,
                                                                       41168, 41448, 54256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 64768, 0, 3,
                                                                       62338, 52672, 62608,
                                                                       42008, 42288, 54616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 65218, 0, 3,
                                                                       62608, 52888, 62878,
                                                                       42288, 42568, 54976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 65668, 0, 3,
                                                                       62878, 53104, 63148,
                                                                       42568, 42848, 55336,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 66118, 0, 3,
                                                                       63418, 53536, 63868,
                                                                       43408, 43828, 55696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 66793, 0, 3,
                                                                       63868, 53896, 64318,
                                                                       43828, 44248, 56236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 67468, 0, 3,
                                                                       64768, 54616, 65218,
                                                                       45088, 45508, 56776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 68143, 0, 3,
                                                                       65218, 54976, 65668,
                                                                       45508, 45928, 57316,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 68818, 0, 3,
                                                                       66118, 55696, 66793,
                                                                       46768, 47356, 57856,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 69763, 0, 3,
                                                                       67468, 56776, 68143,
                                                                       48532, 49120, 58612,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 70708, 69763, 945, ncols);

                    simdfunc::contract_primitives(buffer, 71653, 68818, 945, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 72598, 70708, 21, 1, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 72598, 17, nmax);

        simdtrf::transform_l_inner(buffer, 72598, 71653, 21, 1, nmax);

        simdtrf::transform_h_outer(values + 187 * nvalues + n * npairs, nvalues, buffer, 72598,
                                   17, nmax);
    }

    for (size_t m = 0; m < 374; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
