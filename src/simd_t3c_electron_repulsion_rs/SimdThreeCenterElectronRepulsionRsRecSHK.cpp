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


#include "SimdThreeCenterElectronRepulsionRsRecSHK.hpp"

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
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_shk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_shk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 49959, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 330 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 49959, 48132, 1512, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12}, ncols,
                                                            fj, i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 19, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 77, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 80, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 83, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 86, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 89, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 92, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 95, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 98, 0, 3, 7, 8,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 104, 0, 3, 8, 9,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 110, 0, 3, 9, 10,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 116, 0, 3, 10, 11,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 122, 0, 3, 11, 12,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 128, 0, 3, 12, 13,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 134, 0, 3, 13, 14,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 140, 0, 3, 14, 15,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 146, 0, 3, 15, 16,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 152, 0, 3, 16, 17,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 158, 0, 3, 20, 21,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 164, 0, 3, 21, 22,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 170, 0, 3, 22, 23,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 176, 0, 3, 23, 24,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 182, 0, 3, 24, 25,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 188, 0, 3, 25, 26,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 194, 0, 3, 26, 27,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 200, 0, 3, 27, 28,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 206, 0, 3, 28, 29,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 212, 0, 3, 29, 30,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 32, 35,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 35, 38,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 38, 41,
                                                                       110, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 41, 44,
                                                                       116, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 44, 47,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 47, 50,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 50, 53,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 53, 56,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 56, 59,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 65, 68,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 68, 71,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 71, 74,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 74, 77,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 77, 80,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 80, 83,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 83, 86,
                                                                       194, 200, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 86, 89,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 89, 92,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 98,
                                                                       104, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 413, 0, 3, 104,
                                                                       110, 228, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 110,
                                                                       116, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 443, 0, 3, 116,
                                                                       122, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 122,
                                                                       128, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 473, 0, 3, 128,
                                                                       134, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 134,
                                                                       140, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 503, 0, 3, 140,
                                                                       146, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 158,
                                                                       164, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 533, 0, 3, 164,
                                                                       170, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 170,
                                                                       176, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 563, 0, 3, 176,
                                                                       182, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 578, 0, 3, 182,
                                                                       188, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 593, 0, 3, 188,
                                                                       194, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 608, 0, 3, 194,
                                                                       200, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 623, 0, 3, 200,
                                                                       206, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 638, 0, 3, 218,
                                                                       228, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 659, 0, 3, 228,
                                                                       238, 413, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 680, 0, 3, 238,
                                                                       248, 428, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 701, 0, 3, 248,
                                                                       258, 443, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 722, 0, 3, 258,
                                                                       268, 458, 473, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 743, 0, 3, 268,
                                                                       278, 473, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 764, 0, 3, 278,
                                                                       288, 488, 503, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 785, 0, 3, 308,
                                                                       318, 518, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 806, 0, 3, 318,
                                                                       328, 533, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 827, 0, 3, 328,
                                                                       338, 548, 563, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 848, 0, 3, 338,
                                                                       348, 563, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 869, 0, 3, 348,
                                                                       358, 578, 593, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 890, 0, 3, 358,
                                                                       368, 593, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 911, 0, 3, 368,
                                                                       378, 608, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 932, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 935, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 938, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 941, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 944, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 947, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 950, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 953, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 956, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 959, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 962, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 965, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 968, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 971, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 974, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 977, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 980, 3, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 983, 3, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 986, 3, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 989, 3, 27, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 992, 3, 28, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 995, 3, 29, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 998, 3, 30, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1001, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1004, 3, 9, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1013, 3, 10, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1022, 3, 11, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1031, 3, 12, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1040, 3, 13, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1049, 3, 14, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1058, 3, 15, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1067, 3, 16, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1076, 3, 17, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1085, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1094, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1103, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1112, 3, 25, 80,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1121, 3, 26, 83,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1130, 3, 27, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1139, 3, 28, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1148, 3, 29, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1157, 3, 30, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1166, 3, 32, 98,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1184, 3, 35, 104,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1202, 3, 38, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1220, 3, 41, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1238, 3, 44, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1256, 3, 47, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1274, 3, 50, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1292, 3, 53, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1310, 3, 56, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1328, 3, 59, 152,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1346, 3, 65, 158,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1364, 3, 68, 164,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1382, 3, 71, 170,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1400, 3, 74, 176,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1418, 3, 77, 182,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1436, 3, 80, 188,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1454, 3, 83, 194,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1472, 3, 86, 200,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1490, 3, 89, 206,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1508, 3, 92, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1526, 3, 98, 218,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1556, 3, 104, 228,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1586, 3, 110, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1616, 3, 116, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1646, 3, 122, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1676, 3, 128, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1706, 3, 134, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1736, 3, 140, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1766, 3, 146, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1796, 3, 158, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1826, 3, 164, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1856, 3, 170, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1886, 3, 176, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1916, 3, 182, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1946, 3, 188, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1976, 3, 194, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2006, 3, 200, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2036, 3, 206, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2066, 3, 218, 398,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2111, 3, 228, 413,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2156, 3, 238, 428,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2201, 3, 248, 443,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2246, 3, 258, 458,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2291, 3, 268, 473,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2336, 3, 278, 488,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2381, 3, 288, 503,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2426, 3, 308, 518,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2471, 3, 318, 533,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2516, 3, 328, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2561, 3, 338, 563,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2606, 3, 348, 578,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2651, 3, 358, 593,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2696, 3, 368, 608,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2741, 3, 378, 623,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2786, 3, 398, 638,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2849, 3, 413, 659,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2912, 3, 428, 680,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2975, 3, 443, 701,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3038, 3, 458, 722,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3101, 3, 473, 743,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3164, 3, 488, 764,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3227, 3, 518, 785,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3290, 3, 533, 806,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3353, 3, 548, 827,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3416, 3, 563, 848,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3479, 3, 578, 869,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3542, 3, 593, 890,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3605, 3, 608, 911,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3668, 3, 7, 8,
                                                                       938, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3674, 3, 8, 9,
                                                                       941, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3680, 3, 9, 10,
                                                                       944, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3686, 3, 10, 11,
                                                                       947, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3692, 3, 11, 12,
                                                                       950, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3698, 3, 12, 13,
                                                                       953, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3704, 3, 13, 14,
                                                                       956, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3710, 3, 14, 15,
                                                                       959, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3716, 3, 15, 16,
                                                                       962, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3722, 3, 16, 17,
                                                                       965, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3728, 3, 20, 21,
                                                                       974, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3734, 3, 21, 22,
                                                                       977, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3740, 3, 22, 23,
                                                                       980, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3746, 3, 23, 24,
                                                                       983, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3752, 3, 24, 25,
                                                                       986, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3758, 3, 25, 26,
                                                                       989, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3764, 3, 26, 27,
                                                                       992, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3770, 3, 27, 28,
                                                                       995, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3776, 3, 28, 29,
                                                                       998, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3782, 3, 29, 30,
                                                                       1001, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3788, 0, 3, 3668,
                                                                       938, 3674, 1004, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3806, 0, 3, 3674,
                                                                       941, 3680, 1013, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3824, 0, 3, 3680,
                                                                       944, 3686, 1022, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3842, 0, 3, 3686,
                                                                       947, 3692, 1031, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3860, 0, 3, 3692,
                                                                       950, 3698, 1040, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3878, 0, 3, 3698,
                                                                       953, 3704, 1049, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3896, 0, 3, 3704,
                                                                       956, 3710, 1058, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3914, 0, 3, 3710,
                                                                       959, 3716, 1067, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3932, 0, 3, 3716,
                                                                       962, 3722, 1076, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3950, 0, 3, 3728,
                                                                       974, 3734, 1085, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 3734,
                                                                       977, 3740, 1094, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3986, 0, 3, 3740,
                                                                       980, 3746, 1103, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4004, 0, 3, 3746,
                                                                       983, 3752, 1112, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4022, 0, 3, 3752,
                                                                       986, 3758, 1121, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4040, 0, 3, 3758,
                                                                       989, 3764, 1130, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4058, 0, 3, 3764,
                                                                       992, 3770, 1139, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4076, 0, 3, 3770,
                                                                       995, 3776, 1148, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4094, 0, 3, 3776,
                                                                       998, 3782, 1157, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4112, 0, 3, 3788,
                                                                       1004, 3806, 98, 104, 1202,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4148, 0, 3, 3806,
                                                                       1013, 3824, 104, 110,
                                                                       1220, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4184, 0, 3, 3824,
                                                                       1022, 3842, 110, 116,
                                                                       1238, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4220, 0, 3, 3842,
                                                                       1031, 3860, 116, 122,
                                                                       1256, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4256, 0, 3, 3860,
                                                                       1040, 3878, 122, 128,
                                                                       1274, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4292, 0, 3, 3878,
                                                                       1049, 3896, 128, 134,
                                                                       1292, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4328, 0, 3, 3896,
                                                                       1058, 3914, 134, 140,
                                                                       1310, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4364, 0, 3, 3914,
                                                                       1067, 3932, 140, 146,
                                                                       1328, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4400, 0, 3, 3950,
                                                                       1085, 3968, 158, 164,
                                                                       1382, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4436, 0, 3, 3968,
                                                                       1094, 3986, 164, 170,
                                                                       1400, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4472, 0, 3, 3986,
                                                                       1103, 4004, 170, 176,
                                                                       1418, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4508, 0, 3, 4004,
                                                                       1112, 4022, 176, 182,
                                                                       1436, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4544, 0, 3, 4022,
                                                                       1121, 4040, 182, 188,
                                                                       1454, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4580, 0, 3, 4040,
                                                                       1130, 4058, 188, 194,
                                                                       1472, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4616, 0, 3, 4058,
                                                                       1139, 4076, 194, 200,
                                                                       1490, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4652, 0, 3, 4076,
                                                                       1148, 4094, 200, 206,
                                                                       1508, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4688, 0, 3, 4112,
                                                                       1202, 4148, 218, 228,
                                                                       1586, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4748, 0, 3, 4148,
                                                                       1220, 4184, 228, 238,
                                                                       1616, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4808, 0, 3, 4184,
                                                                       1238, 4220, 238, 248,
                                                                       1646, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4868, 0, 3, 4220,
                                                                       1256, 4256, 248, 258,
                                                                       1676, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4928, 0, 3, 4256,
                                                                       1274, 4292, 258, 268,
                                                                       1706, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4988, 0, 3, 4292,
                                                                       1292, 4328, 268, 278,
                                                                       1736, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5048, 0, 3, 4328,
                                                                       1310, 4364, 278, 288,
                                                                       1766, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5108, 0, 3, 4400,
                                                                       1382, 4436, 308, 318,
                                                                       1856, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5168, 0, 3, 4436,
                                                                       1400, 4472, 318, 328,
                                                                       1886, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5228, 0, 3, 4472,
                                                                       1418, 4508, 328, 338,
                                                                       1916, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 4508,
                                                                       1436, 4544, 338, 348,
                                                                       1946, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5348, 0, 3, 4544,
                                                                       1454, 4580, 348, 358,
                                                                       1976, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5408, 0, 3, 4580,
                                                                       1472, 4616, 358, 368,
                                                                       2006, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5468, 0, 3, 4616,
                                                                       1490, 4652, 368, 378,
                                                                       2036, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5528, 0, 3, 4688,
                                                                       1586, 4748, 398, 413,
                                                                       2156, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5618, 0, 3, 4748,
                                                                       1616, 4808, 413, 428,
                                                                       2201, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5708, 0, 3, 4808,
                                                                       1646, 4868, 428, 443,
                                                                       2246, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5798, 0, 3, 4868,
                                                                       1676, 4928, 443, 458,
                                                                       2291, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5888, 0, 3, 4928,
                                                                       1706, 4988, 458, 473,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5978, 0, 3, 4988,
                                                                       1736, 5048, 473, 488,
                                                                       2381, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6068, 0, 3, 5108,
                                                                       1856, 5168, 518, 533,
                                                                       2516, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6158, 0, 3, 5168,
                                                                       1886, 5228, 533, 548,
                                                                       2561, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 5228,
                                                                       1916, 5288, 548, 563,
                                                                       2606, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6338, 0, 3, 5288,
                                                                       1946, 5348, 563, 578,
                                                                       2651, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 5348,
                                                                       1976, 5408, 578, 593,
                                                                       2696, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6518, 0, 3, 5408,
                                                                       2006, 5468, 593, 608,
                                                                       2741, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 5528,
                                                                       2156, 5618, 638, 659,
                                                                       2912, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6734, 0, 3, 5618,
                                                                       2201, 5708, 659, 680,
                                                                       2975, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6860, 0, 3, 5708,
                                                                       2246, 5798, 680, 701,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6986, 0, 3, 5798,
                                                                       2291, 5888, 701, 722,
                                                                       3101, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7112, 0, 3, 5888,
                                                                       2336, 5978, 722, 743,
                                                                       3164, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7238, 0, 3, 6068,
                                                                       2516, 6158, 785, 806,
                                                                       3353, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7364, 0, 3, 6158,
                                                                       2561, 6248, 806, 827,
                                                                       3416, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7490, 0, 3, 6248,
                                                                       2606, 6338, 827, 848,
                                                                       3479, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7616, 0, 3, 6338,
                                                                       2651, 6428, 848, 869,
                                                                       3542, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7742, 0, 3, 6428,
                                                                       2696, 6518, 869, 890,
                                                                       3605, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7868, 3, 932, 935,
                                                                       3668, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7878, 3, 935, 938,
                                                                       3674, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7888, 3, 938, 941,
                                                                       3680, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7898, 3, 941, 944,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7908, 3, 944, 947,
                                                                       3692, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7918, 3, 947, 950,
                                                                       3698, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7928, 3, 950, 953,
                                                                       3704, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7938, 3, 953, 956,
                                                                       3710, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7948, 3, 956, 959,
                                                                       3716, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7958, 3, 959, 962,
                                                                       3722, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7968, 3, 968, 971,
                                                                       3728, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7978, 3, 971, 974,
                                                                       3734, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7988, 3, 974, 977,
                                                                       3740, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7998, 3, 977, 980,
                                                                       3746, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8008, 3, 980, 983,
                                                                       3752, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8018, 3, 983, 986,
                                                                       3758, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8028, 3, 986, 989,
                                                                       3764, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8038, 3, 989, 992,
                                                                       3770, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8048, 3, 992, 995,
                                                                       3776, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8058, 3, 995, 998,
                                                                       3782, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8068, 0, 3, 7868,
                                                                       3668, 7878, 3788, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8098, 0, 3, 7878,
                                                                       3674, 7888, 3806, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8128, 0, 3, 7888,
                                                                       3680, 7898, 3824, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8158, 0, 3, 7898,
                                                                       3686, 7908, 3842, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8188, 0, 3, 7908,
                                                                       3692, 7918, 3860, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8218, 0, 3, 7918,
                                                                       3698, 7928, 3878, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8248, 0, 3, 7928,
                                                                       3704, 7938, 3896, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8278, 0, 3, 7938,
                                                                       3710, 7948, 3914, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8308, 0, 3, 7948,
                                                                       3716, 7958, 3932, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8338, 0, 3, 7968,
                                                                       3728, 7978, 3950, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8368, 0, 3, 7978,
                                                                       3734, 7988, 3968, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8398, 0, 3, 7988,
                                                                       3740, 7998, 3986, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8428, 0, 3, 7998,
                                                                       3746, 8008, 4004, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8458, 0, 3, 8008,
                                                                       3752, 8018, 4022, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8488, 0, 3, 8018,
                                                                       3758, 8028, 4040, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8518, 0, 3, 8028,
                                                                       3764, 8038, 4058, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8548, 0, 3, 8038,
                                                                       3770, 8048, 4076, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8578, 0, 3, 8048,
                                                                       3776, 8058, 4094, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8608, 0, 3, 8068,
                                                                       3788, 8098, 1166, 1184,
                                                                       4112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8668, 0, 3, 8098,
                                                                       3806, 8128, 1184, 1202,
                                                                       4148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8728, 0, 3, 8128,
                                                                       3824, 8158, 1202, 1220,
                                                                       4184, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8788, 0, 3, 8158,
                                                                       3842, 8188, 1220, 1238,
                                                                       4220, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8848, 0, 3, 8188,
                                                                       3860, 8218, 1238, 1256,
                                                                       4256, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8908, 0, 3, 8218,
                                                                       3878, 8248, 1256, 1274,
                                                                       4292, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8968, 0, 3, 8248,
                                                                       3896, 8278, 1274, 1292,
                                                                       4328, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9028, 0, 3, 8278,
                                                                       3914, 8308, 1292, 1310,
                                                                       4364, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9088, 0, 3, 8338,
                                                                       3950, 8368, 1346, 1364,
                                                                       4400, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9148, 0, 3, 8368,
                                                                       3968, 8398, 1364, 1382,
                                                                       4436, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9208, 0, 3, 8398,
                                                                       3986, 8428, 1382, 1400,
                                                                       4472, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9268, 0, 3, 8428,
                                                                       4004, 8458, 1400, 1418,
                                                                       4508, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9328, 0, 3, 8458,
                                                                       4022, 8488, 1418, 1436,
                                                                       4544, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9388, 0, 3, 8488,
                                                                       4040, 8518, 1436, 1454,
                                                                       4580, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9448, 0, 3, 8518,
                                                                       4058, 8548, 1454, 1472,
                                                                       4616, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9508, 0, 3, 8548,
                                                                       4076, 8578, 1472, 1490,
                                                                       4652, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9568, 0, 3, 8608,
                                                                       4112, 8668, 1526, 1556,
                                                                       4688, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9668, 0, 3, 8668,
                                                                       4148, 8728, 1556, 1586,
                                                                       4748, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9768, 0, 3, 8728,
                                                                       4184, 8788, 1586, 1616,
                                                                       4808, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9868, 0, 3, 8788,
                                                                       4220, 8848, 1616, 1646,
                                                                       4868, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9968, 0, 3, 8848,
                                                                       4256, 8908, 1646, 1676,
                                                                       4928, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10068, 0, 3, 8908,
                                                                       4292, 8968, 1676, 1706,
                                                                       4988, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10168, 0, 3, 8968,
                                                                       4328, 9028, 1706, 1736,
                                                                       5048, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10268, 0, 3, 9088,
                                                                       4400, 9148, 1796, 1826,
                                                                       5108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10368, 0, 3, 9148,
                                                                       4436, 9208, 1826, 1856,
                                                                       5168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10468, 0, 3, 9208,
                                                                       4472, 9268, 1856, 1886,
                                                                       5228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10568, 0, 3, 9268,
                                                                       4508, 9328, 1886, 1916,
                                                                       5288, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10668, 0, 3, 9328,
                                                                       4544, 9388, 1916, 1946,
                                                                       5348, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10768, 0, 3, 9388,
                                                                       4580, 9448, 1946, 1976,
                                                                       5408, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10868, 0, 3, 9448,
                                                                       4616, 9508, 1976, 2006,
                                                                       5468, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10968, 0, 3, 9568,
                                                                       4688, 9668, 2066, 2111,
                                                                       5528, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11118, 0, 3, 9668,
                                                                       4748, 9768, 2111, 2156,
                                                                       5618, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11268, 0, 3, 9768,
                                                                       4808, 9868, 2156, 2201,
                                                                       5708, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11418, 0, 3, 9868,
                                                                       4868, 9968, 2201, 2246,
                                                                       5798, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11568, 0, 3, 9968,
                                                                       4928, 10068, 2246, 2291,
                                                                       5888, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11718, 0, 3,
                                                                       10068, 4988, 10168, 2291,
                                                                       2336, 5978, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11868, 0, 3,
                                                                       10268, 5108, 10368, 2426,
                                                                       2471, 6068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12018, 0, 3,
                                                                       10368, 5168, 10468, 2471,
                                                                       2516, 6158, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12168, 0, 3,
                                                                       10468, 5228, 10568, 2516,
                                                                       2561, 6248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12318, 0, 3,
                                                                       10568, 5288, 10668, 2561,
                                                                       2606, 6338, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12468, 0, 3,
                                                                       10668, 5348, 10768, 2606,
                                                                       2651, 6428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12618, 0, 3,
                                                                       10768, 5408, 10868, 2651,
                                                                       2696, 6518, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12768, 0, 3,
                                                                       10968, 5528, 11118, 2786,
                                                                       2849, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12978, 0, 3,
                                                                       11118, 5618, 11268, 2849,
                                                                       2912, 6734, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13188, 0, 3,
                                                                       11268, 5708, 11418, 2912,
                                                                       2975, 6860, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13398, 0, 3,
                                                                       11418, 5798, 11568, 2975,
                                                                       3038, 6986, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13608, 0, 3,
                                                                       11568, 5888, 11718, 3038,
                                                                       3101, 7112, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13818, 0, 3,
                                                                       11868, 6068, 12018, 3227,
                                                                       3290, 7238, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14028, 0, 3,
                                                                       12018, 6158, 12168, 3290,
                                                                       3353, 7364, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14238, 0, 3,
                                                                       12168, 6248, 12318, 3353,
                                                                       3416, 7490, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14448, 0, 3,
                                                                       12318, 6338, 12468, 3416,
                                                                       3479, 7616, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14658, 0, 3,
                                                                       12468, 6428, 12618, 3479,
                                                                       3542, 7742, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14868, 3, 3668,
                                                                       3674, 7888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14883, 3, 3674,
                                                                       3680, 7898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14898, 3, 3680,
                                                                       3686, 7908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14913, 3, 3686,
                                                                       3692, 7918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14928, 3, 3692,
                                                                       3698, 7928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14943, 3, 3698,
                                                                       3704, 7938, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14958, 3, 3704,
                                                                       3710, 7948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14973, 3, 3710,
                                                                       3716, 7958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14988, 3, 3728,
                                                                       3734, 7988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15003, 3, 3734,
                                                                       3740, 7998, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15018, 3, 3740,
                                                                       3746, 8008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15033, 3, 3746,
                                                                       3752, 8018, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15048, 3, 3752,
                                                                       3758, 8028, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15063, 3, 3758,
                                                                       3764, 8038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15078, 3, 3764,
                                                                       3770, 8048, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15093, 3, 3770,
                                                                       3776, 8058, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15108, 0, 3,
                                                                       14868, 7888, 14883, 3788,
                                                                       3806, 8128, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15153, 0, 3,
                                                                       14883, 7898, 14898, 3806,
                                                                       3824, 8158, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15198, 0, 3,
                                                                       14898, 7908, 14913, 3824,
                                                                       3842, 8188, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15243, 0, 3,
                                                                       14913, 7918, 14928, 3842,
                                                                       3860, 8218, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15288, 0, 3,
                                                                       14928, 7928, 14943, 3860,
                                                                       3878, 8248, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15333, 0, 3,
                                                                       14943, 7938, 14958, 3878,
                                                                       3896, 8278, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15378, 0, 3,
                                                                       14958, 7948, 14973, 3896,
                                                                       3914, 8308, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15423, 0, 3,
                                                                       14988, 7988, 15003, 3950,
                                                                       3968, 8398, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15468, 0, 3,
                                                                       15003, 7998, 15018, 3968,
                                                                       3986, 8428, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15513, 0, 3,
                                                                       15018, 8008, 15033, 3986,
                                                                       4004, 8458, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15558, 0, 3,
                                                                       15033, 8018, 15048, 4004,
                                                                       4022, 8488, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15603, 0, 3,
                                                                       15048, 8028, 15063, 4022,
                                                                       4040, 8518, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15648, 0, 3,
                                                                       15063, 8038, 15078, 4040,
                                                                       4058, 8548, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15693, 0, 3,
                                                                       15078, 8048, 15093, 4058,
                                                                       4076, 8578, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 15738, 0, 3,
                                                                       15108, 8128, 15153, 4112,
                                                                       4148, 8728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 15828, 0, 3,
                                                                       15153, 8158, 15198, 4148,
                                                                       4184, 8788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 15918, 0, 3,
                                                                       15198, 8188, 15243, 4184,
                                                                       4220, 8848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 16008, 0, 3,
                                                                       15243, 8218, 15288, 4220,
                                                                       4256, 8908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 16098, 0, 3,
                                                                       15288, 8248, 15333, 4256,
                                                                       4292, 8968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 16188, 0, 3,
                                                                       15333, 8278, 15378, 4292,
                                                                       4328, 9028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 16278, 0, 3,
                                                                       15423, 8398, 15468, 4400,
                                                                       4436, 9208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 16368, 0, 3,
                                                                       15468, 8428, 15513, 4436,
                                                                       4472, 9268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 16458, 0, 3,
                                                                       15513, 8458, 15558, 4472,
                                                                       4508, 9328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 16548, 0, 3,
                                                                       15558, 8488, 15603, 4508,
                                                                       4544, 9388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 16638, 0, 3,
                                                                       15603, 8518, 15648, 4544,
                                                                       4580, 9448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 16728, 0, 3,
                                                                       15648, 8548, 15693, 4580,
                                                                       4616, 9508, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 16818, 0, 3,
                                                                       15738, 8728, 15828, 4688,
                                                                       4748, 9768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 16968, 0, 3,
                                                                       15828, 8788, 15918, 4748,
                                                                       4808, 9868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 17118, 0, 3,
                                                                       15918, 8848, 16008, 4808,
                                                                       4868, 9968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 17268, 0, 3,
                                                                       16008, 8908, 16098, 4868,
                                                                       4928, 10068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 17418, 0, 3,
                                                                       16098, 8968, 16188, 4928,
                                                                       4988, 10168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 17568, 0, 3,
                                                                       16278, 9208, 16368, 5108,
                                                                       5168, 10468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 17718, 0, 3,
                                                                       16368, 9268, 16458, 5168,
                                                                       5228, 10568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 17868, 0, 3,
                                                                       16458, 9328, 16548, 5228,
                                                                       5288, 10668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18018, 0, 3,
                                                                       16548, 9388, 16638, 5288,
                                                                       5348, 10768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18168, 0, 3,
                                                                       16638, 9448, 16728, 5348,
                                                                       5408, 10868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 18318, 0, 3,
                                                                       16818, 9768, 16968, 5528,
                                                                       5618, 11268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 18543, 0, 3,
                                                                       16968, 9868, 17118, 5618,
                                                                       5708, 11418, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 18768, 0, 3,
                                                                       17118, 9968, 17268, 5708,
                                                                       5798, 11568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 18993, 0, 3,
                                                                       17268, 10068, 17418, 5798,
                                                                       5888, 11718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 19218, 0, 3,
                                                                       17568, 10468, 17718, 6068,
                                                                       6158, 12168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 19443, 0, 3,
                                                                       17718, 10568, 17868, 6158,
                                                                       6248, 12318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 19668, 0, 3,
                                                                       17868, 10668, 18018, 6248,
                                                                       6338, 12468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 19893, 0, 3,
                                                                       18018, 10768, 18168, 6338,
                                                                       6428, 12618, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 20118, 0, 3,
                                                                       18318, 11268, 18543, 6608,
                                                                       6734, 13188, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 20433, 0, 3,
                                                                       18543, 11418, 18768, 6734,
                                                                       6860, 13398, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 20748, 0, 3,
                                                                       18768, 11568, 18993, 6860,
                                                                       6986, 13608, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 21063, 0, 3,
                                                                       19218, 12168, 19443, 7238,
                                                                       7364, 14238, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 21378, 0, 3,
                                                                       19443, 12318, 19668, 7364,
                                                                       7490, 14448, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 21693, 0, 3,
                                                                       19668, 12468, 19893, 7490,
                                                                       7616, 14658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22008, 3, 7868,
                                                                       7878, 14868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22029, 3, 7878,
                                                                       7888, 14883, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22050, 3, 7888,
                                                                       7898, 14898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22071, 3, 7898,
                                                                       7908, 14913, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22092, 3, 7908,
                                                                       7918, 14928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22113, 3, 7918,
                                                                       7928, 14943, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22134, 3, 7928,
                                                                       7938, 14958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22155, 3, 7938,
                                                                       7948, 14973, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22176, 3, 7968,
                                                                       7978, 14988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22197, 3, 7978,
                                                                       7988, 15003, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22218, 3, 7988,
                                                                       7998, 15018, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22239, 3, 7998,
                                                                       8008, 15033, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22260, 3, 8008,
                                                                       8018, 15048, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22281, 3, 8018,
                                                                       8028, 15063, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22302, 3, 8028,
                                                                       8038, 15078, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22323, 3, 8038,
                                                                       8048, 15093, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22344, 0, 3,
                                                                       22008, 14868, 22029, 8068,
                                                                       8098, 15108, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22407, 0, 3,
                                                                       22029, 14883, 22050, 8098,
                                                                       8128, 15153, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22470, 0, 3,
                                                                       22050, 14898, 22071, 8128,
                                                                       8158, 15198, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22533, 0, 3,
                                                                       22071, 14913, 22092, 8158,
                                                                       8188, 15243, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22596, 0, 3,
                                                                       22092, 14928, 22113, 8188,
                                                                       8218, 15288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22659, 0, 3,
                                                                       22113, 14943, 22134, 8218,
                                                                       8248, 15333, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22722, 0, 3,
                                                                       22134, 14958, 22155, 8248,
                                                                       8278, 15378, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22785, 0, 3,
                                                                       22176, 14988, 22197, 8338,
                                                                       8368, 15423, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22848, 0, 3,
                                                                       22197, 15003, 22218, 8368,
                                                                       8398, 15468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22911, 0, 3,
                                                                       22218, 15018, 22239, 8398,
                                                                       8428, 15513, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 22974, 0, 3,
                                                                       22239, 15033, 22260, 8428,
                                                                       8458, 15558, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 23037, 0, 3,
                                                                       22260, 15048, 22281, 8458,
                                                                       8488, 15603, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 23100, 0, 3,
                                                                       22281, 15063, 22302, 8488,
                                                                       8518, 15648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 23163, 0, 3,
                                                                       22302, 15078, 22323, 8518,
                                                                       8548, 15693, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 23226, 0, 3,
                                                                       22344, 15108, 22407, 8608,
                                                                       8668, 15738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 23352, 0, 3,
                                                                       22407, 15153, 22470, 8668,
                                                                       8728, 15828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 23478, 0, 3,
                                                                       22470, 15198, 22533, 8728,
                                                                       8788, 15918, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 23604, 0, 3,
                                                                       22533, 15243, 22596, 8788,
                                                                       8848, 16008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 23730, 0, 3,
                                                                       22596, 15288, 22659, 8848,
                                                                       8908, 16098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 23856, 0, 3,
                                                                       22659, 15333, 22722, 8908,
                                                                       8968, 16188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 23982, 0, 3,
                                                                       22785, 15423, 22848, 9088,
                                                                       9148, 16278, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 24108, 0, 3,
                                                                       22848, 15468, 22911, 9148,
                                                                       9208, 16368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 24234, 0, 3,
                                                                       22911, 15513, 22974, 9208,
                                                                       9268, 16458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 24360, 0, 3,
                                                                       22974, 15558, 23037, 9268,
                                                                       9328, 16548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 24486, 0, 3,
                                                                       23037, 15603, 23100, 9328,
                                                                       9388, 16638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 24612, 0, 3,
                                                                       23100, 15648, 23163, 9388,
                                                                       9448, 16728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 24738, 0, 3,
                                                                       23226, 15738, 23352, 9568,
                                                                       9668, 16818, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 24948, 0, 3,
                                                                       23352, 15828, 23478, 9668,
                                                                       9768, 16968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 25158, 0, 3,
                                                                       23478, 15918, 23604, 9768,
                                                                       9868, 17118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 25368, 0, 3,
                                                                       23604, 16008, 23730, 9868,
                                                                       9968, 17268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 25578, 0, 3,
                                                                       23730, 16098, 23856, 9968,
                                                                       10068, 17418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 25788, 0, 3,
                                                                       23982, 16278, 24108,
                                                                       10268, 10368, 17568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 25998, 0, 3,
                                                                       24108, 16368, 24234,
                                                                       10368, 10468, 17718,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 26208, 0, 3,
                                                                       24234, 16458, 24360,
                                                                       10468, 10568, 17868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 26418, 0, 3,
                                                                       24360, 16548, 24486,
                                                                       10568, 10668, 18018,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 26628, 0, 3,
                                                                       24486, 16638, 24612,
                                                                       10668, 10768, 18168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 26838, 0, 3,
                                                                       24738, 16818, 24948,
                                                                       10968, 11118, 18318,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 27153, 0, 3,
                                                                       24948, 16968, 25158,
                                                                       11118, 11268, 18543,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 27468, 0, 3,
                                                                       25158, 17118, 25368,
                                                                       11268, 11418, 18768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 27783, 0, 3,
                                                                       25368, 17268, 25578,
                                                                       11418, 11568, 18993,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 28098, 0, 3,
                                                                       25788, 17568, 25998,
                                                                       11868, 12018, 19218,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 28413, 0, 3,
                                                                       25998, 17718, 26208,
                                                                       12018, 12168, 19443,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 28728, 0, 3,
                                                                       26208, 17868, 26418,
                                                                       12168, 12318, 19668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 29043, 0, 3,
                                                                       26418, 18018, 26628,
                                                                       12318, 12468, 19893,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 29358, 0, 3,
                                                                       26838, 18318, 27153,
                                                                       12768, 12978, 20118,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 29799, 0, 3,
                                                                       27153, 18543, 27468,
                                                                       12978, 13188, 20433,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 30240, 0, 3,
                                                                       27468, 18768, 27783,
                                                                       13188, 13398, 20748,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 30681, 0, 3,
                                                                       28098, 19218, 28413,
                                                                       13818, 14028, 21063,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 31122, 0, 3,
                                                                       28413, 19443, 28728,
                                                                       14028, 14238, 21378,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 31563, 0, 3,
                                                                       28728, 19668, 29043,
                                                                       14238, 14448, 21693,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32004, 3, 14868,
                                                                       14883, 22050, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32032, 3, 14883,
                                                                       14898, 22071, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32060, 3, 14898,
                                                                       14913, 22092, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32088, 3, 14913,
                                                                       14928, 22113, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32116, 3, 14928,
                                                                       14943, 22134, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32144, 3, 14943,
                                                                       14958, 22155, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32172, 3, 14988,
                                                                       15003, 22218, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32200, 3, 15003,
                                                                       15018, 22239, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32228, 3, 15018,
                                                                       15033, 22260, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32256, 3, 15033,
                                                                       15048, 22281, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32284, 3, 15048,
                                                                       15063, 22302, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32312, 3, 15063,
                                                                       15078, 22323, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 32340, 0, 3,
                                                                       32004, 22050, 32032,
                                                                       15108, 15153, 22470,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 32424, 0, 3,
                                                                       32032, 22071, 32060,
                                                                       15153, 15198, 22533,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 32508, 0, 3,
                                                                       32060, 22092, 32088,
                                                                       15198, 15243, 22596,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 32592, 0, 3,
                                                                       32088, 22113, 32116,
                                                                       15243, 15288, 22659,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 32676, 0, 3,
                                                                       32116, 22134, 32144,
                                                                       15288, 15333, 22722,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 32760, 0, 3,
                                                                       32172, 22218, 32200,
                                                                       15423, 15468, 22911,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 32844, 0, 3,
                                                                       32200, 22239, 32228,
                                                                       15468, 15513, 22974,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 32928, 0, 3,
                                                                       32228, 22260, 32256,
                                                                       15513, 15558, 23037,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33012, 0, 3,
                                                                       32256, 22281, 32284,
                                                                       15558, 15603, 23100,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33096, 0, 3,
                                                                       32284, 22302, 32312,
                                                                       15603, 15648, 23163,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 33180, 0, 3,
                                                                       32340, 22470, 32424,
                                                                       15738, 15828, 23478,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 33348, 0, 3,
                                                                       32424, 22533, 32508,
                                                                       15828, 15918, 23604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 33516, 0, 3,
                                                                       32508, 22596, 32592,
                                                                       15918, 16008, 23730,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 33684, 0, 3,
                                                                       32592, 22659, 32676,
                                                                       16008, 16098, 23856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 33852, 0, 3,
                                                                       32760, 22911, 32844,
                                                                       16278, 16368, 24234,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 34020, 0, 3,
                                                                       32844, 22974, 32928,
                                                                       16368, 16458, 24360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 34188, 0, 3,
                                                                       32928, 23037, 33012,
                                                                       16458, 16548, 24486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 34356, 0, 3,
                                                                       33012, 23100, 33096,
                                                                       16548, 16638, 24612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 34524, 0, 3,
                                                                       33180, 23478, 33348,
                                                                       16818, 16968, 25158,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 34804, 0, 3,
                                                                       33348, 23604, 33516,
                                                                       16968, 17118, 25368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 35084, 0, 3,
                                                                       33516, 23730, 33684,
                                                                       17118, 17268, 25578,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 35364, 0, 3,
                                                                       33852, 24234, 34020,
                                                                       17568, 17718, 26208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 35644, 0, 3,
                                                                       34020, 24360, 34188,
                                                                       17718, 17868, 26418,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 35924, 0, 3,
                                                                       34188, 24486, 34356,
                                                                       17868, 18018, 26628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 36204, 0, 3,
                                                                       34524, 25158, 34804,
                                                                       18318, 18543, 27468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 36624, 0, 3,
                                                                       34804, 25368, 35084,
                                                                       18543, 18768, 27783,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 37044, 0, 3,
                                                                       35364, 26208, 35644,
                                                                       19218, 19443, 28728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 37464, 0, 3,
                                                                       35644, 26418, 35924,
                                                                       19443, 19668, 29043,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 37884, 0, 3,
                                                                       36204, 27468, 36624,
                                                                       20118, 20433, 30240,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 38472, 0, 3,
                                                                       37044, 28728, 37464,
                                                                       21063, 21378, 31563,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39060, 3, 22008,
                                                                       22029, 32004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39096, 3, 22029,
                                                                       22050, 32032, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39132, 3, 22050,
                                                                       22071, 32060, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39168, 3, 22071,
                                                                       22092, 32088, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39204, 3, 22092,
                                                                       22113, 32116, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39240, 3, 22113,
                                                                       22134, 32144, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39276, 3, 22176,
                                                                       22197, 32172, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39312, 3, 22197,
                                                                       22218, 32200, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39348, 3, 22218,
                                                                       22239, 32228, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39384, 3, 22239,
                                                                       22260, 32256, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39420, 3, 22260,
                                                                       22281, 32284, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39456, 3, 22281,
                                                                       22302, 32312, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 39492, 0, 3,
                                                                       39060, 32004, 39096,
                                                                       22344, 22407, 32340,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 39600, 0, 3,
                                                                       39096, 32032, 39132,
                                                                       22407, 22470, 32424,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 39708, 0, 3,
                                                                       39132, 32060, 39168,
                                                                       22470, 22533, 32508,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 39816, 0, 3,
                                                                       39168, 32088, 39204,
                                                                       22533, 22596, 32592,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 39924, 0, 3,
                                                                       39204, 32116, 39240,
                                                                       22596, 22659, 32676,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 40032, 0, 3,
                                                                       39276, 32172, 39312,
                                                                       22785, 22848, 32760,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 40140, 0, 3,
                                                                       39312, 32200, 39348,
                                                                       22848, 22911, 32844,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 40248, 0, 3,
                                                                       39348, 32228, 39384,
                                                                       22911, 22974, 32928,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 40356, 0, 3,
                                                                       39384, 32256, 39420,
                                                                       22974, 23037, 33012,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 40464, 0, 3,
                                                                       39420, 32284, 39456,
                                                                       23037, 23100, 33096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 40572, 0, 3,
                                                                       39492, 32340, 39600,
                                                                       23226, 23352, 33180,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 40788, 0, 3,
                                                                       39600, 32424, 39708,
                                                                       23352, 23478, 33348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 41004, 0, 3,
                                                                       39708, 32508, 39816,
                                                                       23478, 23604, 33516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 41220, 0, 3,
                                                                       39816, 32592, 39924,
                                                                       23604, 23730, 33684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 41436, 0, 3,
                                                                       40032, 32760, 40140,
                                                                       23982, 24108, 33852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 41652, 0, 3,
                                                                       40140, 32844, 40248,
                                                                       24108, 24234, 34020,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 41868, 0, 3,
                                                                       40248, 32928, 40356,
                                                                       24234, 24360, 34188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 42084, 0, 3,
                                                                       40356, 33012, 40464,
                                                                       24360, 24486, 34356,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 42300, 0, 3,
                                                                       40572, 33180, 40788,
                                                                       24738, 24948, 34524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 42660, 0, 3,
                                                                       40788, 33348, 41004,
                                                                       24948, 25158, 34804,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 43020, 0, 3,
                                                                       41004, 33516, 41220,
                                                                       25158, 25368, 35084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 43380, 0, 3,
                                                                       41436, 33852, 41652,
                                                                       25788, 25998, 35364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 43740, 0, 3,
                                                                       41652, 34020, 41868,
                                                                       25998, 26208, 35644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 44100, 0, 3,
                                                                       41868, 34188, 42084,
                                                                       26208, 26418, 35924,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 44460, 0, 3,
                                                                       42300, 34524, 42660,
                                                                       26838, 27153, 36204,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 45000, 0, 3,
                                                                       42660, 34804, 43020,
                                                                       27153, 27468, 36624,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 45540, 0, 3,
                                                                       43380, 35364, 43740,
                                                                       28098, 28413, 37044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 46080, 0, 3,
                                                                       43740, 35644, 44100,
                                                                       28413, 28728, 37464,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 46620, 0, 3,
                                                                       44460, 36204, 45000,
                                                                       29358, 29799, 37884,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 47376, 0, 3,
                                                                       45540, 37044, 46080,
                                                                       30681, 31122, 38472,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 48132, 47376, 756, ncols);

                    simdfunc::contract_primitives(buffer, 48888, 46620, 756, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 49644, 48132, 21, 1, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 49644, 15, nmax);

        simdtrf::transform_k_inner(buffer, 49644, 48888, 21, 1, nmax);

        simdtrf::transform_h_outer(values + 165 * nvalues + n * npairs, nvalues, buffer, 49644,
                                   15, nmax);
    }

    for (size_t m = 0; m < 330; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
