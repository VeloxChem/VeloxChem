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


#include "SimdThreeCenterElectronRepulsionRsRecPIG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_pig_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_pig_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 33323, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 702 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 33323, 28388, 2748, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 11,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 19, 3, 11,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 932, 0, 3, 398,
                                                                       413, 638, 659, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 960, 0, 3, 413,
                                                                       428, 659, 680, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 988, 0, 3, 428,
                                                                       443, 680, 701, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1016, 0, 3, 443,
                                                                       458, 701, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 458,
                                                                       473, 722, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 473,
                                                                       488, 743, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 518,
                                                                       533, 785, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 533,
                                                                       548, 806, 827, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 548,
                                                                       563, 827, 848, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 563,
                                                                       578, 848, 869, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 578,
                                                                       593, 869, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 593,
                                                                       608, 890, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 638,
                                                                       659, 932, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1304, 0, 3, 659,
                                                                       680, 960, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1340, 0, 3, 680,
                                                                       701, 988, 1016, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1376, 0, 3, 701,
                                                                       722, 1016, 1044, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1412, 0, 3, 722,
                                                                       743, 1044, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1448, 0, 3, 785,
                                                                       806, 1100, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1484, 0, 3, 806,
                                                                       827, 1128, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 827,
                                                                       848, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1556, 0, 3, 848,
                                                                       869, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1592, 0, 3, 869,
                                                                       890, 1212, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1628, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1631, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1634, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1637, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1640, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1643, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1646, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1649, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1652, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1655, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1658, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1661, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1664, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1667, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1670, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1673, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1676, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1679, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1682, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1685, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1688, 3, 9, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1697, 3, 10, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1706, 3, 11, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1715, 3, 12, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1724, 3, 13, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1733, 3, 14, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1742, 3, 15, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1751, 3, 16, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1760, 3, 17, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1769, 3, 22, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1778, 3, 23, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1787, 3, 24, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1796, 3, 25, 80,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1805, 3, 26, 83,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1814, 3, 27, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1823, 3, 28, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1832, 3, 29, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1841, 3, 30, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1850, 3, 38, 110,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1868, 3, 41, 116,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1886, 3, 44, 122,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1904, 3, 47, 128,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1922, 3, 50, 134,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1940, 3, 53, 140,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1958, 3, 56, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1976, 3, 59, 152,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1994, 3, 71, 170,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2012, 3, 74, 176,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2030, 3, 77, 182,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2048, 3, 80, 188,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2066, 3, 83, 194,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2084, 3, 86, 200,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2102, 3, 89, 206,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2120, 3, 92, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2138, 3, 110, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2168, 3, 116, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2198, 3, 122, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2228, 3, 128, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2258, 3, 134, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2288, 3, 140, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2318, 3, 146, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2348, 3, 170, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2378, 3, 176, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2408, 3, 182, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2438, 3, 188, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2468, 3, 194, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2498, 3, 200, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2528, 3, 206, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2558, 3, 238, 428,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2603, 3, 248, 443,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2648, 3, 258, 458,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2693, 3, 268, 473,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2738, 3, 278, 488,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2783, 3, 288, 503,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2828, 3, 328, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2873, 3, 338, 563,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2918, 3, 348, 578,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2963, 3, 358, 593,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3008, 3, 368, 608,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3053, 3, 378, 623,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3098, 3, 428, 680,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3161, 3, 443, 701,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3224, 3, 458, 722,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3287, 3, 473, 743,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3350, 3, 488, 764,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3413, 3, 548, 827,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3476, 3, 563, 848,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3539, 3, 578, 869,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3602, 3, 593, 890,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3665, 3, 608, 911,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3728, 3, 680, 988,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3812, 3, 701,
                                                                       1016, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3896, 3, 722,
                                                                       1044, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3980, 3, 743,
                                                                       1072, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4064, 3, 827,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4148, 3, 848,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4232, 3, 869,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4316, 3, 890,
                                                                       1240, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4400, 3, 988,
                                                                       1340, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4508, 3, 1016,
                                                                       1376, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4616, 3, 1044,
                                                                       1412, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4724, 3, 1156,
                                                                       1520, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4832, 3, 1184,
                                                                       1556, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4940, 3, 1212,
                                                                       1592, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5048, 3, 7, 8,
                                                                       1628, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5054, 3, 8, 9,
                                                                       1631, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5060, 3, 9, 10,
                                                                       1634, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5066, 3, 10, 11,
                                                                       1637, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5072, 3, 11, 12,
                                                                       1640, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5078, 3, 12, 13,
                                                                       1643, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5084, 3, 13, 14,
                                                                       1646, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5090, 3, 14, 15,
                                                                       1649, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5096, 3, 15, 16,
                                                                       1652, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5102, 3, 16, 17,
                                                                       1655, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5108, 3, 20, 21,
                                                                       1658, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5114, 3, 21, 22,
                                                                       1661, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5120, 3, 22, 23,
                                                                       1664, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5126, 3, 23, 24,
                                                                       1667, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5132, 3, 24, 25,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5138, 3, 25, 26,
                                                                       1673, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5144, 3, 26, 27,
                                                                       1676, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5150, 3, 27, 28,
                                                                       1679, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5156, 3, 28, 29,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5162, 3, 29, 30,
                                                                       1685, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5168, 0, 3, 5048,
                                                                       1628, 5054, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5186, 0, 3, 5054,
                                                                       1631, 5060, 1697, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5204, 0, 3, 5060,
                                                                       1634, 5066, 1706, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5222, 0, 3, 5066,
                                                                       1637, 5072, 1715, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5240, 0, 3, 5072,
                                                                       1640, 5078, 1724, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5258, 0, 3, 5078,
                                                                       1643, 5084, 1733, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5276, 0, 3, 5084,
                                                                       1646, 5090, 1742, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5294, 0, 3, 5090,
                                                                       1649, 5096, 1751, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5312, 0, 3, 5096,
                                                                       1652, 5102, 1760, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5330, 0, 3, 5108,
                                                                       1658, 5114, 1769, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5348, 0, 3, 5114,
                                                                       1661, 5120, 1778, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5366, 0, 3, 5120,
                                                                       1664, 5126, 1787, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5384, 0, 3, 5126,
                                                                       1667, 5132, 1796, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5402, 0, 3, 5132,
                                                                       1670, 5138, 1805, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5420, 0, 3, 5138,
                                                                       1673, 5144, 1814, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5438, 0, 3, 5144,
                                                                       1676, 5150, 1823, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5456, 0, 3, 5150,
                                                                       1679, 5156, 1832, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5474, 0, 3, 5156,
                                                                       1682, 5162, 1841, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5492, 0, 3, 5168,
                                                                       1688, 5186, 98, 104, 1850,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5528, 0, 3, 5186,
                                                                       1697, 5204, 104, 110,
                                                                       1868, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5564, 0, 3, 5204,
                                                                       1706, 5222, 110, 116,
                                                                       1886, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5600, 0, 3, 5222,
                                                                       1715, 5240, 116, 122,
                                                                       1904, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5636, 0, 3, 5240,
                                                                       1724, 5258, 122, 128,
                                                                       1922, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5672, 0, 3, 5258,
                                                                       1733, 5276, 128, 134,
                                                                       1940, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5708, 0, 3, 5276,
                                                                       1742, 5294, 134, 140,
                                                                       1958, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5744, 0, 3, 5294,
                                                                       1751, 5312, 140, 146,
                                                                       1976, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5780, 0, 3, 5330,
                                                                       1769, 5348, 158, 164,
                                                                       1994, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5816, 0, 3, 5348,
                                                                       1778, 5366, 164, 170,
                                                                       2012, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5852, 0, 3, 5366,
                                                                       1787, 5384, 170, 176,
                                                                       2030, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5888, 0, 3, 5384,
                                                                       1796, 5402, 176, 182,
                                                                       2048, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5924, 0, 3, 5402,
                                                                       1805, 5420, 182, 188,
                                                                       2066, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5960, 0, 3, 5420,
                                                                       1814, 5438, 188, 194,
                                                                       2084, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5996, 0, 3, 5438,
                                                                       1823, 5456, 194, 200,
                                                                       2102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6032, 0, 3, 5456,
                                                                       1832, 5474, 200, 206,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6068, 0, 3, 5492,
                                                                       1850, 5528, 218, 228,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6128, 0, 3, 5528,
                                                                       1868, 5564, 228, 238,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6188, 0, 3, 5564,
                                                                       1886, 5600, 238, 248,
                                                                       2198, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 5600,
                                                                       1904, 5636, 248, 258,
                                                                       2228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6308, 0, 3, 5636,
                                                                       1922, 5672, 258, 268,
                                                                       2258, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6368, 0, 3, 5672,
                                                                       1940, 5708, 268, 278,
                                                                       2288, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 5708,
                                                                       1958, 5744, 278, 288,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6488, 0, 3, 5780,
                                                                       1994, 5816, 308, 318,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6548, 0, 3, 5816,
                                                                       2012, 5852, 318, 328,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 5852,
                                                                       2030, 5888, 328, 338,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6668, 0, 3, 5888,
                                                                       2048, 5924, 338, 348,
                                                                       2438, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6728, 0, 3, 5924,
                                                                       2066, 5960, 348, 358,
                                                                       2468, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 5960,
                                                                       2084, 5996, 358, 368,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6848, 0, 3, 5996,
                                                                       2102, 6032, 368, 378,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6908, 0, 3, 6068,
                                                                       2138, 6128, 398, 413,
                                                                       2558, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6998, 0, 3, 6128,
                                                                       2168, 6188, 413, 428,
                                                                       2603, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7088, 0, 3, 6188,
                                                                       2198, 6248, 428, 443,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7178, 0, 3, 6248,
                                                                       2228, 6308, 443, 458,
                                                                       2693, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7268, 0, 3, 6308,
                                                                       2258, 6368, 458, 473,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7358, 0, 3, 6368,
                                                                       2288, 6428, 473, 488,
                                                                       2783, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7448, 0, 3, 6488,
                                                                       2348, 6548, 518, 533,
                                                                       2828, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7538, 0, 3, 6548,
                                                                       2378, 6608, 533, 548,
                                                                       2873, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7628, 0, 3, 6608,
                                                                       2408, 6668, 548, 563,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7718, 0, 3, 6668,
                                                                       2438, 6728, 563, 578,
                                                                       2963, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7808, 0, 3, 6728,
                                                                       2468, 6788, 578, 593,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7898, 0, 3, 6788,
                                                                       2498, 6848, 593, 608,
                                                                       3053, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7988, 0, 3, 6908,
                                                                       2558, 6998, 638, 659,
                                                                       3098, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8114, 0, 3, 6998,
                                                                       2603, 7088, 659, 680,
                                                                       3161, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8240, 0, 3, 7088,
                                                                       2648, 7178, 680, 701,
                                                                       3224, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8366, 0, 3, 7178,
                                                                       2693, 7268, 701, 722,
                                                                       3287, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8492, 0, 3, 7268,
                                                                       2738, 7358, 722, 743,
                                                                       3350, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8618, 0, 3, 7448,
                                                                       2828, 7538, 785, 806,
                                                                       3413, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8744, 0, 3, 7538,
                                                                       2873, 7628, 806, 827,
                                                                       3476, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8870, 0, 3, 7628,
                                                                       2918, 7718, 827, 848,
                                                                       3539, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8996, 0, 3, 7718,
                                                                       2963, 7808, 848, 869,
                                                                       3602, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9122, 0, 3, 7808,
                                                                       3008, 7898, 869, 890,
                                                                       3665, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9248, 0, 3, 7988,
                                                                       3098, 8114, 932, 960,
                                                                       3728, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9416, 0, 3, 8114,
                                                                       3161, 8240, 960, 988,
                                                                       3812, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9584, 0, 3, 8240,
                                                                       3224, 8366, 988, 1016,
                                                                       3896, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9752, 0, 3, 8366,
                                                                       3287, 8492, 1016, 1044,
                                                                       3980, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9920, 0, 3, 8618,
                                                                       3413, 8744, 1100, 1128,
                                                                       4064, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10088, 0, 3, 8744,
                                                                       3476, 8870, 1128, 1156,
                                                                       4148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10256, 0, 3, 8870,
                                                                       3539, 8996, 1156, 1184,
                                                                       4232, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10424, 0, 3, 8996,
                                                                       3602, 9122, 1184, 1212,
                                                                       4316, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10592, 0, 3, 9248,
                                                                       3728, 9416, 1268, 1304,
                                                                       4400, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10808, 0, 3, 9416,
                                                                       3812, 9584, 1304, 1340,
                                                                       4508, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11024, 0, 3, 9584,
                                                                       3896, 9752, 1340, 1376,
                                                                       4616, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11240, 0, 3, 9920,
                                                                       4064, 10088, 1448, 1484,
                                                                       4724, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11456, 0, 3,
                                                                       10088, 4148, 10256, 1484,
                                                                       1520, 4832, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11672, 0, 3,
                                                                       10256, 4232, 10424, 1520,
                                                                       1556, 4940, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11888, 3, 1628,
                                                                       1631, 5060, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11898, 3, 1631,
                                                                       1634, 5066, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11908, 3, 1634,
                                                                       1637, 5072, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11918, 3, 1637,
                                                                       1640, 5078, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11928, 3, 1640,
                                                                       1643, 5084, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11938, 3, 1643,
                                                                       1646, 5090, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11948, 3, 1646,
                                                                       1649, 5096, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11958, 3, 1649,
                                                                       1652, 5102, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11968, 3, 1658,
                                                                       1661, 5120, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11978, 3, 1661,
                                                                       1664, 5126, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11988, 3, 1664,
                                                                       1667, 5132, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11998, 3, 1667,
                                                                       1670, 5138, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12008, 3, 1670,
                                                                       1673, 5144, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12018, 3, 1673,
                                                                       1676, 5150, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12028, 3, 1676,
                                                                       1679, 5156, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12038, 3, 1679,
                                                                       1682, 5162, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12048, 0, 3,
                                                                       11888, 5060, 11898, 5204,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12078, 0, 3,
                                                                       11898, 5066, 11908, 5222,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12108, 0, 3,
                                                                       11908, 5072, 11918, 5240,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12138, 0, 3,
                                                                       11918, 5078, 11928, 5258,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12168, 0, 3,
                                                                       11928, 5084, 11938, 5276,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12198, 0, 3,
                                                                       11938, 5090, 11948, 5294,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12228, 0, 3,
                                                                       11948, 5096, 11958, 5312,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12258, 0, 3,
                                                                       11968, 5120, 11978, 5366,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12288, 0, 3,
                                                                       11978, 5126, 11988, 5384,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12318, 0, 3,
                                                                       11988, 5132, 11998, 5402,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12348, 0, 3,
                                                                       11998, 5138, 12008, 5420,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12378, 0, 3,
                                                                       12008, 5144, 12018, 5438,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12408, 0, 3,
                                                                       12018, 5150, 12028, 5456,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12438, 0, 3,
                                                                       12028, 5156, 12038, 5474,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12468, 0, 3,
                                                                       12048, 5204, 12078, 1850,
                                                                       1868, 5564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12528, 0, 3,
                                                                       12078, 5222, 12108, 1868,
                                                                       1886, 5600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12588, 0, 3,
                                                                       12108, 5240, 12138, 1886,
                                                                       1904, 5636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12648, 0, 3,
                                                                       12138, 5258, 12168, 1904,
                                                                       1922, 5672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12708, 0, 3,
                                                                       12168, 5276, 12198, 1922,
                                                                       1940, 5708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12768, 0, 3,
                                                                       12198, 5294, 12228, 1940,
                                                                       1958, 5744, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12828, 0, 3,
                                                                       12258, 5366, 12288, 1994,
                                                                       2012, 5852, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12888, 0, 3,
                                                                       12288, 5384, 12318, 2012,
                                                                       2030, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12948, 0, 3,
                                                                       12318, 5402, 12348, 2030,
                                                                       2048, 5924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13008, 0, 3,
                                                                       12348, 5420, 12378, 2048,
                                                                       2066, 5960, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13068, 0, 3,
                                                                       12378, 5438, 12408, 2066,
                                                                       2084, 5996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13128, 0, 3,
                                                                       12408, 5456, 12438, 2084,
                                                                       2102, 6032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13188, 0, 3,
                                                                       12468, 5564, 12528, 2138,
                                                                       2168, 6188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13288, 0, 3,
                                                                       12528, 5600, 12588, 2168,
                                                                       2198, 6248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13388, 0, 3,
                                                                       12588, 5636, 12648, 2198,
                                                                       2228, 6308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13488, 0, 3,
                                                                       12648, 5672, 12708, 2228,
                                                                       2258, 6368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13588, 0, 3,
                                                                       12708, 5708, 12768, 2258,
                                                                       2288, 6428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13688, 0, 3,
                                                                       12828, 5852, 12888, 2348,
                                                                       2378, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13788, 0, 3,
                                                                       12888, 5888, 12948, 2378,
                                                                       2408, 6668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13888, 0, 3,
                                                                       12948, 5924, 13008, 2408,
                                                                       2438, 6728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13988, 0, 3,
                                                                       13008, 5960, 13068, 2438,
                                                                       2468, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14088, 0, 3,
                                                                       13068, 5996, 13128, 2468,
                                                                       2498, 6848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14188, 0, 3,
                                                                       13188, 6188, 13288, 2558,
                                                                       2603, 7088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14338, 0, 3,
                                                                       13288, 6248, 13388, 2603,
                                                                       2648, 7178, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14488, 0, 3,
                                                                       13388, 6308, 13488, 2648,
                                                                       2693, 7268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14638, 0, 3,
                                                                       13488, 6368, 13588, 2693,
                                                                       2738, 7358, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14788, 0, 3,
                                                                       13688, 6608, 13788, 2828,
                                                                       2873, 7628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14938, 0, 3,
                                                                       13788, 6668, 13888, 2873,
                                                                       2918, 7718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15088, 0, 3,
                                                                       13888, 6728, 13988, 2918,
                                                                       2963, 7808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15238, 0, 3,
                                                                       13988, 6788, 14088, 2963,
                                                                       3008, 7898, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15388, 0, 3,
                                                                       14188, 7088, 14338, 3098,
                                                                       3161, 8240, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15598, 0, 3,
                                                                       14338, 7178, 14488, 3161,
                                                                       3224, 8366, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15808, 0, 3,
                                                                       14488, 7268, 14638, 3224,
                                                                       3287, 8492, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16018, 0, 3,
                                                                       14788, 7628, 14938, 3413,
                                                                       3476, 8870, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16228, 0, 3,
                                                                       14938, 7718, 15088, 3476,
                                                                       3539, 8996, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16438, 0, 3,
                                                                       15088, 7808, 15238, 3539,
                                                                       3602, 9122, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16648, 0, 3,
                                                                       15388, 8240, 15598, 3728,
                                                                       3812, 9584, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16928, 0, 3,
                                                                       15598, 8366, 15808, 3812,
                                                                       3896, 9752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 17208, 0, 3,
                                                                       16018, 8870, 16228, 4064,
                                                                       4148, 10256, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 17488, 0, 3,
                                                                       16228, 8996, 16438, 4148,
                                                                       4232, 10424, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 17768, 0, 3,
                                                                       16648, 9584, 16928, 4400,
                                                                       4508, 11024, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 18128, 0, 3,
                                                                       17208, 10256, 17488, 4724,
                                                                       4832, 11672, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18488, 3, 5048,
                                                                       5054, 11888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18503, 3, 5054,
                                                                       5060, 11898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18518, 3, 5060,
                                                                       5066, 11908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18533, 3, 5066,
                                                                       5072, 11918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18548, 3, 5072,
                                                                       5078, 11928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18563, 3, 5078,
                                                                       5084, 11938, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18578, 3, 5084,
                                                                       5090, 11948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18593, 3, 5090,
                                                                       5096, 11958, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18608, 3, 5108,
                                                                       5114, 11968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18623, 3, 5114,
                                                                       5120, 11978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18638, 3, 5120,
                                                                       5126, 11988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18653, 3, 5126,
                                                                       5132, 11998, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18668, 3, 5132,
                                                                       5138, 12008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18683, 3, 5138,
                                                                       5144, 12018, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18698, 3, 5144,
                                                                       5150, 12028, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18713, 3, 5150,
                                                                       5156, 12038, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18728, 0, 3,
                                                                       18488, 11888, 18503, 5168,
                                                                       5186, 12048, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18773, 0, 3,
                                                                       18503, 11898, 18518, 5186,
                                                                       5204, 12078, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18818, 0, 3,
                                                                       18518, 11908, 18533, 5204,
                                                                       5222, 12108, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18863, 0, 3,
                                                                       18533, 11918, 18548, 5222,
                                                                       5240, 12138, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18908, 0, 3,
                                                                       18548, 11928, 18563, 5240,
                                                                       5258, 12168, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18953, 0, 3,
                                                                       18563, 11938, 18578, 5258,
                                                                       5276, 12198, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18998, 0, 3,
                                                                       18578, 11948, 18593, 5276,
                                                                       5294, 12228, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19043, 0, 3,
                                                                       18608, 11968, 18623, 5330,
                                                                       5348, 12258, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19088, 0, 3,
                                                                       18623, 11978, 18638, 5348,
                                                                       5366, 12288, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19133, 0, 3,
                                                                       18638, 11988, 18653, 5366,
                                                                       5384, 12318, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19178, 0, 3,
                                                                       18653, 11998, 18668, 5384,
                                                                       5402, 12348, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19223, 0, 3,
                                                                       18668, 12008, 18683, 5402,
                                                                       5420, 12378, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19268, 0, 3,
                                                                       18683, 12018, 18698, 5420,
                                                                       5438, 12408, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19313, 0, 3,
                                                                       18698, 12028, 18713, 5438,
                                                                       5456, 12438, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19358, 0, 3,
                                                                       18728, 12048, 18773, 5492,
                                                                       5528, 12468, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19448, 0, 3,
                                                                       18773, 12078, 18818, 5528,
                                                                       5564, 12528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19538, 0, 3,
                                                                       18818, 12108, 18863, 5564,
                                                                       5600, 12588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19628, 0, 3,
                                                                       18863, 12138, 18908, 5600,
                                                                       5636, 12648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19718, 0, 3,
                                                                       18908, 12168, 18953, 5636,
                                                                       5672, 12708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19808, 0, 3,
                                                                       18953, 12198, 18998, 5672,
                                                                       5708, 12768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19898, 0, 3,
                                                                       19043, 12258, 19088, 5780,
                                                                       5816, 12828, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19988, 0, 3,
                                                                       19088, 12288, 19133, 5816,
                                                                       5852, 12888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20078, 0, 3,
                                                                       19133, 12318, 19178, 5852,
                                                                       5888, 12948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20168, 0, 3,
                                                                       19178, 12348, 19223, 5888,
                                                                       5924, 13008, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20258, 0, 3,
                                                                       19223, 12378, 19268, 5924,
                                                                       5960, 13068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20348, 0, 3,
                                                                       19268, 12408, 19313, 5960,
                                                                       5996, 13128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20438, 0, 3,
                                                                       19358, 12468, 19448, 6068,
                                                                       6128, 13188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20588, 0, 3,
                                                                       19448, 12528, 19538, 6128,
                                                                       6188, 13288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20738, 0, 3,
                                                                       19538, 12588, 19628, 6188,
                                                                       6248, 13388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20888, 0, 3,
                                                                       19628, 12648, 19718, 6248,
                                                                       6308, 13488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21038, 0, 3,
                                                                       19718, 12708, 19808, 6308,
                                                                       6368, 13588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21188, 0, 3,
                                                                       19898, 12828, 19988, 6488,
                                                                       6548, 13688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21338, 0, 3,
                                                                       19988, 12888, 20078, 6548,
                                                                       6608, 13788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21488, 0, 3,
                                                                       20078, 12948, 20168, 6608,
                                                                       6668, 13888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21638, 0, 3,
                                                                       20168, 13008, 20258, 6668,
                                                                       6728, 13988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21788, 0, 3,
                                                                       20258, 13068, 20348, 6728,
                                                                       6788, 14088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21938, 0, 3,
                                                                       20438, 13188, 20588, 6908,
                                                                       6998, 14188, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22163, 0, 3,
                                                                       20588, 13288, 20738, 6998,
                                                                       7088, 14338, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22388, 0, 3,
                                                                       20738, 13388, 20888, 7088,
                                                                       7178, 14488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22613, 0, 3,
                                                                       20888, 13488, 21038, 7178,
                                                                       7268, 14638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22838, 0, 3,
                                                                       21188, 13688, 21338, 7448,
                                                                       7538, 14788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23063, 0, 3,
                                                                       21338, 13788, 21488, 7538,
                                                                       7628, 14938, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23288, 0, 3,
                                                                       21488, 13888, 21638, 7628,
                                                                       7718, 15088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23513, 0, 3,
                                                                       21638, 13988, 21788, 7718,
                                                                       7808, 15238, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 23738, 0, 3,
                                                                       21938, 14188, 22163, 7988,
                                                                       8114, 15388, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24053, 0, 3,
                                                                       22163, 14338, 22388, 8114,
                                                                       8240, 15598, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24368, 0, 3,
                                                                       22388, 14488, 22613, 8240,
                                                                       8366, 15808, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24683, 0, 3,
                                                                       22838, 14788, 23063, 8618,
                                                                       8744, 16018, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24998, 0, 3,
                                                                       23063, 14938, 23288, 8744,
                                                                       8870, 16228, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25313, 0, 3,
                                                                       23288, 15088, 23513, 8870,
                                                                       8996, 16438, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 25628, 0, 3,
                                                                       23738, 15388, 24053, 9248,
                                                                       9416, 16648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 26048, 0, 3,
                                                                       24053, 15598, 24368, 9416,
                                                                       9584, 16928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 26468, 0, 3,
                                                                       24683, 16018, 24998, 9920,
                                                                       10088, 17208, ncols,
                                                                       gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 26888, 0, 3,
                                                                       24998, 16228, 25313,
                                                                       10088, 10256, 17488,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 27308, 0, 3,
                                                                       25628, 16648, 26048,
                                                                       10592, 10808, 17768,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 27848, 0, 3,
                                                                       26468, 17208, 26888,
                                                                       11240, 11456, 18128,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 28388, 25628, 420, ncols);

                    simdfunc::contract_primitives(buffer, 29060, 26468, 420, ncols);

                    simdfunc::contract_primitives(buffer, 29732, 27308, 540, ncols);

                    simdfunc::contract_primitives(buffer, 30596, 27848, 540, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 28808, 28388, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 29480, 29060, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 30272, 29732, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 31136, 30596, 36, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 31460, 28808, 30272, 9, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 32216, 29480, 31136, 9, nmax);

        simdtrf::transform_i_inner(buffer, 32972, 32216, 3, 9, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 32972, 117, nmax);

        simdtrf::transform_i_inner(buffer, 32972, 31460, 3, 9, nmax);

        simdtrf::transform_p_outer(values + 351 * nvalues + n * npairs, nvalues, buffer, 32972,
                                   117, nmax);
    }

    for (size_t m = 0; m < 702; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
