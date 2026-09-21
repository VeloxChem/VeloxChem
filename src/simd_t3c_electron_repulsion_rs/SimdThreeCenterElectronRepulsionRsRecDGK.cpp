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


#include "SimdThreeCenterElectronRepulsionRsRecDGK.hpp"

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
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_dgk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_dgk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 90530, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1350 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 90530, 77252, 6108, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
                                                            ncols, fj, i * nprim_b + j, fq,
                                                            omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 20, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13}, ncols,
                                                        fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 76, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 79, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 82, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 85, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 88, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 91, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 94, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 97, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 100, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 103, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 106, 0, 3, 7, 8,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 112, 0, 3, 8, 9,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 118, 0, 3, 9, 10,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 124, 0, 3, 10, 11,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 130, 0, 3, 11, 12,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 136, 0, 3, 12, 13,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 142, 0, 3, 13, 14,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 148, 0, 3, 14, 15,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 154, 0, 3, 15, 16,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 160, 0, 3, 16, 17,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 166, 0, 3, 17, 18,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 172, 0, 3, 21, 22,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 178, 0, 3, 22, 23,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 184, 0, 3, 23, 24,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 190, 0, 3, 24, 25,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 196, 0, 3, 25, 26,
                                                                       82, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 202, 0, 3, 26, 27,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 208, 0, 3, 27, 28,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 214, 0, 3, 28, 29,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 220, 0, 3, 29, 30,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 226, 0, 3, 30, 31,
                                                                       97, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 232, 0, 3, 31, 32,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 34, 37,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 248, 0, 3, 37, 40,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 40, 43,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 43, 46,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 278, 0, 3, 46, 49,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 49, 52,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 52, 55,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 55, 58,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 58, 61,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 61, 64,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 70, 73,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 73, 76,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 76, 79,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 79, 82,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 82, 85,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 85, 88,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 88, 91,
                                                                       208, 214, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 91, 94,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 94, 97,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 97,
                                                                       100, 226, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 106,
                                                                       112, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 453, 0, 3, 112,
                                                                       118, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 118,
                                                                       124, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 483, 0, 3, 124,
                                                                       130, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 130,
                                                                       136, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 513, 0, 3, 136,
                                                                       142, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 528, 0, 3, 142,
                                                                       148, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 543, 0, 3, 148,
                                                                       154, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 558, 0, 3, 154,
                                                                       160, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 573, 0, 3, 172,
                                                                       178, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 178,
                                                                       184, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 603, 0, 3, 184,
                                                                       190, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 618, 0, 3, 190,
                                                                       196, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 633, 0, 3, 196,
                                                                       202, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 648, 0, 3, 202,
                                                                       208, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 663, 0, 3, 208,
                                                                       214, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 678, 0, 3, 214,
                                                                       220, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 693, 0, 3, 220,
                                                                       226, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 708, 0, 3, 238,
                                                                       248, 438, 453, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 729, 0, 3, 248,
                                                                       258, 453, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 750, 0, 3, 258,
                                                                       268, 468, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 771, 0, 3, 268,
                                                                       278, 483, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 792, 0, 3, 278,
                                                                       288, 498, 513, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 813, 0, 3, 288,
                                                                       298, 513, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 834, 0, 3, 298,
                                                                       308, 528, 543, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 855, 0, 3, 308,
                                                                       318, 543, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 876, 0, 3, 338,
                                                                       348, 573, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 897, 0, 3, 348,
                                                                       358, 588, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 918, 0, 3, 358,
                                                                       368, 603, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 939, 0, 3, 368,
                                                                       378, 618, 633, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 960, 0, 3, 378,
                                                                       388, 633, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 981, 0, 3, 388,
                                                                       398, 648, 663, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 398,
                                                                       408, 663, 678, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1023, 0, 3, 408,
                                                                       418, 678, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 438,
                                                                       453, 708, 729, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 453,
                                                                       468, 729, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 468,
                                                                       483, 750, 771, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 483,
                                                                       498, 771, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 498,
                                                                       513, 792, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 513,
                                                                       528, 813, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 528,
                                                                       543, 834, 855, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 573,
                                                                       588, 876, 897, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 588,
                                                                       603, 897, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 603,
                                                                       618, 918, 939, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 618,
                                                                       633, 939, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 633,
                                                                       648, 960, 981, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 648,
                                                                       663, 981, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 663,
                                                                       678, 1002, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1436, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1439, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1442, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1445, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1448, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1451, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1454, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1457, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1460, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1463, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1466, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1469, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1472, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1475, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1478, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1481, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1484, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1487, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1490, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1493, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1496, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1499, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1502, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1505, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1508, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1511, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1514, 3, 9, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1523, 3, 10, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1532, 3, 11, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1541, 3, 12, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1550, 3, 13, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1559, 3, 14, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1568, 3, 15, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1577, 3, 16, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1586, 3, 17, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1595, 3, 18, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1604, 3, 23, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1613, 3, 24, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1622, 3, 25, 82,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1631, 3, 26, 85,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1640, 3, 27, 88,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1649, 3, 28, 91,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1658, 3, 29, 94,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1667, 3, 30, 97,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1676, 3, 31, 100,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1685, 3, 32, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1694, 3, 34, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1712, 3, 37, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1730, 3, 40, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1748, 3, 43, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1766, 3, 46, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1784, 3, 49, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1802, 3, 52, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1820, 3, 55, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1838, 3, 58, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1856, 3, 61, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1874, 3, 64, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1892, 3, 70, 172,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1910, 3, 73, 178,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1928, 3, 76, 184,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1946, 3, 79, 190,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1964, 3, 82, 196,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1982, 3, 85, 202,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2000, 3, 88, 208,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2018, 3, 91, 214,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2036, 3, 94, 220,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2054, 3, 97, 226,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2072, 3, 100, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2090, 3, 106, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2120, 3, 112, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2150, 3, 118, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2180, 3, 124, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2210, 3, 130, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2240, 3, 136, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2270, 3, 142, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2300, 3, 148, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2330, 3, 154, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2360, 3, 160, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2390, 3, 172, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2420, 3, 178, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2450, 3, 184, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2480, 3, 190, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2510, 3, 196, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2540, 3, 202, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2570, 3, 208, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2600, 3, 214, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2630, 3, 220, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2660, 3, 226, 428,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2690, 3, 238, 438,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2735, 3, 248, 453,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2780, 3, 258, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2825, 3, 268, 483,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2870, 3, 278, 498,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2915, 3, 288, 513,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2960, 3, 298, 528,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3005, 3, 308, 543,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3050, 3, 318, 558,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3095, 3, 338, 573,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3140, 3, 348, 588,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3185, 3, 358, 603,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3230, 3, 368, 618,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3275, 3, 378, 633,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3320, 3, 388, 648,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3365, 3, 398, 663,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3410, 3, 408, 678,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3455, 3, 418, 693,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3500, 3, 438, 708,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3563, 3, 453, 729,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3626, 3, 468, 750,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3689, 3, 483, 771,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3752, 3, 498, 792,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3815, 3, 513, 813,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3878, 3, 528, 834,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3941, 3, 543, 855,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4004, 3, 573, 876,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4067, 3, 588, 897,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4130, 3, 603, 918,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4193, 3, 618, 939,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4256, 3, 633, 960,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4319, 3, 648, 981,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4382, 3, 663,
                                                                       1002, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4445, 3, 678,
                                                                       1023, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4508, 3, 708,
                                                                       1044, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4592, 3, 729,
                                                                       1072, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4676, 3, 750,
                                                                       1100, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4760, 3, 771,
                                                                       1128, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4844, 3, 792,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4928, 3, 813,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5012, 3, 834,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5096, 3, 876,
                                                                       1240, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5180, 3, 897,
                                                                       1268, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5264, 3, 918,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5348, 3, 939,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5432, 3, 960,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5516, 3, 981,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5600, 3, 1002,
                                                                       1408, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5684, 3, 7, 8,
                                                                       1442, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5690, 3, 8, 9,
                                                                       1445, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5696, 3, 9, 10,
                                                                       1448, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5702, 3, 10, 11,
                                                                       1451, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5708, 3, 11, 12,
                                                                       1454, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5714, 3, 12, 13,
                                                                       1457, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5720, 3, 13, 14,
                                                                       1460, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5726, 3, 14, 15,
                                                                       1463, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5732, 3, 15, 16,
                                                                       1466, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5738, 3, 16, 17,
                                                                       1469, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5744, 3, 17, 18,
                                                                       1472, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5750, 3, 21, 22,
                                                                       1481, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5756, 3, 22, 23,
                                                                       1484, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5762, 3, 23, 24,
                                                                       1487, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5768, 3, 24, 25,
                                                                       1490, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5774, 3, 25, 26,
                                                                       1493, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5780, 3, 26, 27,
                                                                       1496, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5786, 3, 27, 28,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5792, 3, 28, 29,
                                                                       1502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5798, 3, 29, 30,
                                                                       1505, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5804, 3, 30, 31,
                                                                       1508, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5810, 3, 31, 32,
                                                                       1511, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5816, 0, 3, 5684,
                                                                       1442, 5690, 1514, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5834, 0, 3, 5690,
                                                                       1445, 5696, 1523, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5852, 0, 3, 5696,
                                                                       1448, 5702, 1532, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5870, 0, 3, 5702,
                                                                       1451, 5708, 1541, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5888, 0, 3, 5708,
                                                                       1454, 5714, 1550, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5906, 0, 3, 5714,
                                                                       1457, 5720, 1559, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5924, 0, 3, 5720,
                                                                       1460, 5726, 1568, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5942, 0, 3, 5726,
                                                                       1463, 5732, 1577, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5960, 0, 3, 5732,
                                                                       1466, 5738, 1586, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5978, 0, 3, 5738,
                                                                       1469, 5744, 1595, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5996, 0, 3, 5750,
                                                                       1481, 5756, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6014, 0, 3, 5756,
                                                                       1484, 5762, 1613, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6032, 0, 3, 5762,
                                                                       1487, 5768, 1622, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6050, 0, 3, 5768,
                                                                       1490, 5774, 1631, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6068, 0, 3, 5774,
                                                                       1493, 5780, 1640, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6086, 0, 3, 5780,
                                                                       1496, 5786, 1649, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6104, 0, 3, 5786,
                                                                       1499, 5792, 1658, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6122, 0, 3, 5792,
                                                                       1502, 5798, 1667, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6140, 0, 3, 5798,
                                                                       1505, 5804, 1676, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6158, 0, 3, 5804,
                                                                       1508, 5810, 1685, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6176, 0, 3, 5816,
                                                                       1514, 5834, 106, 112,
                                                                       1730, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6212, 0, 3, 5834,
                                                                       1523, 5852, 112, 118,
                                                                       1748, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 5852,
                                                                       1532, 5870, 118, 124,
                                                                       1766, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6284, 0, 3, 5870,
                                                                       1541, 5888, 124, 130,
                                                                       1784, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6320, 0, 3, 5888,
                                                                       1550, 5906, 130, 136,
                                                                       1802, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6356, 0, 3, 5906,
                                                                       1559, 5924, 136, 142,
                                                                       1820, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6392, 0, 3, 5924,
                                                                       1568, 5942, 142, 148,
                                                                       1838, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 5942,
                                                                       1577, 5960, 148, 154,
                                                                       1856, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6464, 0, 3, 5960,
                                                                       1586, 5978, 154, 160,
                                                                       1874, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6500, 0, 3, 5996,
                                                                       1604, 6014, 172, 178,
                                                                       1928, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6536, 0, 3, 6014,
                                                                       1613, 6032, 178, 184,
                                                                       1946, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6572, 0, 3, 6032,
                                                                       1622, 6050, 184, 190,
                                                                       1964, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 6050,
                                                                       1631, 6068, 190, 196,
                                                                       1982, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6644, 0, 3, 6068,
                                                                       1640, 6086, 196, 202,
                                                                       2000, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6680, 0, 3, 6086,
                                                                       1649, 6104, 202, 208,
                                                                       2018, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 6104,
                                                                       1658, 6122, 208, 214,
                                                                       2036, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6752, 0, 3, 6122,
                                                                       1667, 6140, 214, 220,
                                                                       2054, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 6140,
                                                                       1676, 6158, 220, 226,
                                                                       2072, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6824, 0, 3, 6176,
                                                                       1730, 6212, 238, 248,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6884, 0, 3, 6212,
                                                                       1748, 6248, 248, 258,
                                                                       2180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6944, 0, 3, 6248,
                                                                       1766, 6284, 258, 268,
                                                                       2210, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7004, 0, 3, 6284,
                                                                       1784, 6320, 268, 278,
                                                                       2240, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7064, 0, 3, 6320,
                                                                       1802, 6356, 278, 288,
                                                                       2270, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7124, 0, 3, 6356,
                                                                       1820, 6392, 288, 298,
                                                                       2300, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7184, 0, 3, 6392,
                                                                       1838, 6428, 298, 308,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7244, 0, 3, 6428,
                                                                       1856, 6464, 308, 318,
                                                                       2360, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7304, 0, 3, 6500,
                                                                       1928, 6536, 338, 348,
                                                                       2450, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7364, 0, 3, 6536,
                                                                       1946, 6572, 348, 358,
                                                                       2480, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7424, 0, 3, 6572,
                                                                       1964, 6608, 358, 368,
                                                                       2510, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7484, 0, 3, 6608,
                                                                       1982, 6644, 368, 378,
                                                                       2540, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7544, 0, 3, 6644,
                                                                       2000, 6680, 378, 388,
                                                                       2570, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7604, 0, 3, 6680,
                                                                       2018, 6716, 388, 398,
                                                                       2600, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7664, 0, 3, 6716,
                                                                       2036, 6752, 398, 408,
                                                                       2630, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7724, 0, 3, 6752,
                                                                       2054, 6788, 408, 418,
                                                                       2660, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7784, 0, 3, 6824,
                                                                       2150, 6884, 438, 453,
                                                                       2780, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7874, 0, 3, 6884,
                                                                       2180, 6944, 453, 468,
                                                                       2825, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7964, 0, 3, 6944,
                                                                       2210, 7004, 468, 483,
                                                                       2870, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8054, 0, 3, 7004,
                                                                       2240, 7064, 483, 498,
                                                                       2915, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8144, 0, 3, 7064,
                                                                       2270, 7124, 498, 513,
                                                                       2960, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8234, 0, 3, 7124,
                                                                       2300, 7184, 513, 528,
                                                                       3005, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8324, 0, 3, 7184,
                                                                       2330, 7244, 528, 543,
                                                                       3050, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8414, 0, 3, 7304,
                                                                       2450, 7364, 573, 588,
                                                                       3185, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8504, 0, 3, 7364,
                                                                       2480, 7424, 588, 603,
                                                                       3230, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8594, 0, 3, 7424,
                                                                       2510, 7484, 603, 618,
                                                                       3275, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8684, 0, 3, 7484,
                                                                       2540, 7544, 618, 633,
                                                                       3320, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8774, 0, 3, 7544,
                                                                       2570, 7604, 633, 648,
                                                                       3365, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8864, 0, 3, 7604,
                                                                       2600, 7664, 648, 663,
                                                                       3410, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8954, 0, 3, 7664,
                                                                       2630, 7724, 663, 678,
                                                                       3455, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9044, 0, 3, 7784,
                                                                       2780, 7874, 708, 729,
                                                                       3626, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9170, 0, 3, 7874,
                                                                       2825, 7964, 729, 750,
                                                                       3689, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9296, 0, 3, 7964,
                                                                       2870, 8054, 750, 771,
                                                                       3752, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9422, 0, 3, 8054,
                                                                       2915, 8144, 771, 792,
                                                                       3815, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9548, 0, 3, 8144,
                                                                       2960, 8234, 792, 813,
                                                                       3878, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9674, 0, 3, 8234,
                                                                       3005, 8324, 813, 834,
                                                                       3941, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9800, 0, 3, 8414,
                                                                       3185, 8504, 876, 897,
                                                                       4130, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9926, 0, 3, 8504,
                                                                       3230, 8594, 897, 918,
                                                                       4193, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10052, 0, 3, 8594,
                                                                       3275, 8684, 918, 939,
                                                                       4256, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10178, 0, 3, 8684,
                                                                       3320, 8774, 939, 960,
                                                                       4319, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10304, 0, 3, 8774,
                                                                       3365, 8864, 960, 981,
                                                                       4382, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10430, 0, 3, 8864,
                                                                       3410, 8954, 981, 1002,
                                                                       4445, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10556, 0, 3, 9044,
                                                                       3626, 9170, 1044, 1072,
                                                                       4676, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10724, 0, 3, 9170,
                                                                       3689, 9296, 1072, 1100,
                                                                       4760, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10892, 0, 3, 9296,
                                                                       3752, 9422, 1100, 1128,
                                                                       4844, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11060, 0, 3, 9422,
                                                                       3815, 9548, 1128, 1156,
                                                                       4928, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11228, 0, 3, 9548,
                                                                       3878, 9674, 1156, 1184,
                                                                       5012, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11396, 0, 3, 9800,
                                                                       4130, 9926, 1240, 1268,
                                                                       5264, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11564, 0, 3, 9926,
                                                                       4193, 10052, 1268, 1296,
                                                                       5348, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11732, 0, 3,
                                                                       10052, 4256, 10178, 1296,
                                                                       1324, 5432, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11900, 0, 3,
                                                                       10178, 4319, 10304, 1324,
                                                                       1352, 5516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12068, 0, 3,
                                                                       10304, 4382, 10430, 1352,
                                                                       1380, 5600, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12236, 3, 1436,
                                                                       1439, 5684, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12246, 3, 1439,
                                                                       1442, 5690, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12256, 3, 1442,
                                                                       1445, 5696, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12266, 3, 1445,
                                                                       1448, 5702, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12276, 3, 1448,
                                                                       1451, 5708, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12286, 3, 1451,
                                                                       1454, 5714, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12296, 3, 1454,
                                                                       1457, 5720, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12306, 3, 1457,
                                                                       1460, 5726, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12316, 3, 1460,
                                                                       1463, 5732, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12326, 3, 1463,
                                                                       1466, 5738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12336, 3, 1466,
                                                                       1469, 5744, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12346, 3, 1475,
                                                                       1478, 5750, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12356, 3, 1478,
                                                                       1481, 5756, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12366, 3, 1481,
                                                                       1484, 5762, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12376, 3, 1484,
                                                                       1487, 5768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12386, 3, 1487,
                                                                       1490, 5774, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12396, 3, 1490,
                                                                       1493, 5780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12406, 3, 1493,
                                                                       1496, 5786, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12416, 3, 1496,
                                                                       1499, 5792, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12426, 3, 1499,
                                                                       1502, 5798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12436, 3, 1502,
                                                                       1505, 5804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12446, 3, 1505,
                                                                       1508, 5810, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12456, 0, 3,
                                                                       12236, 5684, 12246, 5816,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12486, 0, 3,
                                                                       12246, 5690, 12256, 5834,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12516, 0, 3,
                                                                       12256, 5696, 12266, 5852,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12546, 0, 3,
                                                                       12266, 5702, 12276, 5870,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12576, 0, 3,
                                                                       12276, 5708, 12286, 5888,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12606, 0, 3,
                                                                       12286, 5714, 12296, 5906,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12636, 0, 3,
                                                                       12296, 5720, 12306, 5924,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12666, 0, 3,
                                                                       12306, 5726, 12316, 5942,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12696, 0, 3,
                                                                       12316, 5732, 12326, 5960,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12726, 0, 3,
                                                                       12326, 5738, 12336, 5978,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12756, 0, 3,
                                                                       12346, 5750, 12356, 5996,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12786, 0, 3,
                                                                       12356, 5756, 12366, 6014,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12816, 0, 3,
                                                                       12366, 5762, 12376, 6032,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12846, 0, 3,
                                                                       12376, 5768, 12386, 6050,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12876, 0, 3,
                                                                       12386, 5774, 12396, 6068,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12906, 0, 3,
                                                                       12396, 5780, 12406, 6086,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12936, 0, 3,
                                                                       12406, 5786, 12416, 6104,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12966, 0, 3,
                                                                       12416, 5792, 12426, 6122,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12996, 0, 3,
                                                                       12426, 5798, 12436, 6140,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 13026, 0, 3,
                                                                       12436, 5804, 12446, 6158,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13056, 0, 3,
                                                                       12456, 5816, 12486, 1694,
                                                                       1712, 6176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13116, 0, 3,
                                                                       12486, 5834, 12516, 1712,
                                                                       1730, 6212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13176, 0, 3,
                                                                       12516, 5852, 12546, 1730,
                                                                       1748, 6248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13236, 0, 3,
                                                                       12546, 5870, 12576, 1748,
                                                                       1766, 6284, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13296, 0, 3,
                                                                       12576, 5888, 12606, 1766,
                                                                       1784, 6320, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13356, 0, 3,
                                                                       12606, 5906, 12636, 1784,
                                                                       1802, 6356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13416, 0, 3,
                                                                       12636, 5924, 12666, 1802,
                                                                       1820, 6392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13476, 0, 3,
                                                                       12666, 5942, 12696, 1820,
                                                                       1838, 6428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13536, 0, 3,
                                                                       12696, 5960, 12726, 1838,
                                                                       1856, 6464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13596, 0, 3,
                                                                       12756, 5996, 12786, 1892,
                                                                       1910, 6500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13656, 0, 3,
                                                                       12786, 6014, 12816, 1910,
                                                                       1928, 6536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13716, 0, 3,
                                                                       12816, 6032, 12846, 1928,
                                                                       1946, 6572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13776, 0, 3,
                                                                       12846, 6050, 12876, 1946,
                                                                       1964, 6608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13836, 0, 3,
                                                                       12876, 6068, 12906, 1964,
                                                                       1982, 6644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13896, 0, 3,
                                                                       12906, 6086, 12936, 1982,
                                                                       2000, 6680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13956, 0, 3,
                                                                       12936, 6104, 12966, 2000,
                                                                       2018, 6716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14016, 0, 3,
                                                                       12966, 6122, 12996, 2018,
                                                                       2036, 6752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14076, 0, 3,
                                                                       12996, 6140, 13026, 2036,
                                                                       2054, 6788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14136, 0, 3,
                                                                       13056, 6176, 13116, 2090,
                                                                       2120, 6824, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14236, 0, 3,
                                                                       13116, 6212, 13176, 2120,
                                                                       2150, 6884, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14336, 0, 3,
                                                                       13176, 6248, 13236, 2150,
                                                                       2180, 6944, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14436, 0, 3,
                                                                       13236, 6284, 13296, 2180,
                                                                       2210, 7004, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14536, 0, 3,
                                                                       13296, 6320, 13356, 2210,
                                                                       2240, 7064, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14636, 0, 3,
                                                                       13356, 6356, 13416, 2240,
                                                                       2270, 7124, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14736, 0, 3,
                                                                       13416, 6392, 13476, 2270,
                                                                       2300, 7184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14836, 0, 3,
                                                                       13476, 6428, 13536, 2300,
                                                                       2330, 7244, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 14936, 0, 3,
                                                                       13596, 6500, 13656, 2390,
                                                                       2420, 7304, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15036, 0, 3,
                                                                       13656, 6536, 13716, 2420,
                                                                       2450, 7364, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15136, 0, 3,
                                                                       13716, 6572, 13776, 2450,
                                                                       2480, 7424, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15236, 0, 3,
                                                                       13776, 6608, 13836, 2480,
                                                                       2510, 7484, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15336, 0, 3,
                                                                       13836, 6644, 13896, 2510,
                                                                       2540, 7544, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15436, 0, 3,
                                                                       13896, 6680, 13956, 2540,
                                                                       2570, 7604, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15536, 0, 3,
                                                                       13956, 6716, 14016, 2570,
                                                                       2600, 7664, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15636, 0, 3,
                                                                       14016, 6752, 14076, 2600,
                                                                       2630, 7724, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15736, 0, 3,
                                                                       14136, 6824, 14236, 2690,
                                                                       2735, 7784, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 15886, 0, 3,
                                                                       14236, 6884, 14336, 2735,
                                                                       2780, 7874, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16036, 0, 3,
                                                                       14336, 6944, 14436, 2780,
                                                                       2825, 7964, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16186, 0, 3,
                                                                       14436, 7004, 14536, 2825,
                                                                       2870, 8054, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16336, 0, 3,
                                                                       14536, 7064, 14636, 2870,
                                                                       2915, 8144, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16486, 0, 3,
                                                                       14636, 7124, 14736, 2915,
                                                                       2960, 8234, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16636, 0, 3,
                                                                       14736, 7184, 14836, 2960,
                                                                       3005, 8324, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16786, 0, 3,
                                                                       14936, 7304, 15036, 3095,
                                                                       3140, 8414, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16936, 0, 3,
                                                                       15036, 7364, 15136, 3140,
                                                                       3185, 8504, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17086, 0, 3,
                                                                       15136, 7424, 15236, 3185,
                                                                       3230, 8594, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17236, 0, 3,
                                                                       15236, 7484, 15336, 3230,
                                                                       3275, 8684, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17386, 0, 3,
                                                                       15336, 7544, 15436, 3275,
                                                                       3320, 8774, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17536, 0, 3,
                                                                       15436, 7604, 15536, 3320,
                                                                       3365, 8864, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17686, 0, 3,
                                                                       15536, 7664, 15636, 3365,
                                                                       3410, 8954, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 17836, 0, 3,
                                                                       15736, 7784, 15886, 3500,
                                                                       3563, 9044, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 18046, 0, 3,
                                                                       15886, 7874, 16036, 3563,
                                                                       3626, 9170, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 18256, 0, 3,
                                                                       16036, 7964, 16186, 3626,
                                                                       3689, 9296, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 18466, 0, 3,
                                                                       16186, 8054, 16336, 3689,
                                                                       3752, 9422, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 18676, 0, 3,
                                                                       16336, 8144, 16486, 3752,
                                                                       3815, 9548, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 18886, 0, 3,
                                                                       16486, 8234, 16636, 3815,
                                                                       3878, 9674, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19096, 0, 3,
                                                                       16786, 8414, 16936, 4004,
                                                                       4067, 9800, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19306, 0, 3,
                                                                       16936, 8504, 17086, 4067,
                                                                       4130, 9926, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19516, 0, 3,
                                                                       17086, 8594, 17236, 4130,
                                                                       4193, 10052, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19726, 0, 3,
                                                                       17236, 8684, 17386, 4193,
                                                                       4256, 10178, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19936, 0, 3,
                                                                       17386, 8774, 17536, 4256,
                                                                       4319, 10304, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20146, 0, 3,
                                                                       17536, 8864, 17686, 4319,
                                                                       4382, 10430, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 20356, 0, 3,
                                                                       17836, 9044, 18046, 4508,
                                                                       4592, 10556, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 20636, 0, 3,
                                                                       18046, 9170, 18256, 4592,
                                                                       4676, 10724, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 20916, 0, 3,
                                                                       18256, 9296, 18466, 4676,
                                                                       4760, 10892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21196, 0, 3,
                                                                       18466, 9422, 18676, 4760,
                                                                       4844, 11060, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21476, 0, 3,
                                                                       18676, 9548, 18886, 4844,
                                                                       4928, 11228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21756, 0, 3,
                                                                       19096, 9800, 19306, 5096,
                                                                       5180, 11396, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22036, 0, 3,
                                                                       19306, 9926, 19516, 5180,
                                                                       5264, 11564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22316, 0, 3,
                                                                       19516, 10052, 19726, 5264,
                                                                       5348, 11732, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22596, 0, 3,
                                                                       19726, 10178, 19936, 5348,
                                                                       5432, 11900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22876, 0, 3,
                                                                       19936, 10304, 20146, 5432,
                                                                       5516, 12068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23156, 3, 5684,
                                                                       5690, 12256, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23171, 3, 5690,
                                                                       5696, 12266, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23186, 3, 5696,
                                                                       5702, 12276, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23201, 3, 5702,
                                                                       5708, 12286, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23216, 3, 5708,
                                                                       5714, 12296, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23231, 3, 5714,
                                                                       5720, 12306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23246, 3, 5720,
                                                                       5726, 12316, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23261, 3, 5726,
                                                                       5732, 12326, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23276, 3, 5732,
                                                                       5738, 12336, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23291, 3, 5750,
                                                                       5756, 12366, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23306, 3, 5756,
                                                                       5762, 12376, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23321, 3, 5762,
                                                                       5768, 12386, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23336, 3, 5768,
                                                                       5774, 12396, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23351, 3, 5774,
                                                                       5780, 12406, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23366, 3, 5780,
                                                                       5786, 12416, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23381, 3, 5786,
                                                                       5792, 12426, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23396, 3, 5792,
                                                                       5798, 12436, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 23411, 3, 5798,
                                                                       5804, 12446, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23426, 0, 3,
                                                                       23156, 12256, 23171, 5816,
                                                                       5834, 12516, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23471, 0, 3,
                                                                       23171, 12266, 23186, 5834,
                                                                       5852, 12546, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23516, 0, 3,
                                                                       23186, 12276, 23201, 5852,
                                                                       5870, 12576, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23561, 0, 3,
                                                                       23201, 12286, 23216, 5870,
                                                                       5888, 12606, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23606, 0, 3,
                                                                       23216, 12296, 23231, 5888,
                                                                       5906, 12636, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23651, 0, 3,
                                                                       23231, 12306, 23246, 5906,
                                                                       5924, 12666, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23696, 0, 3,
                                                                       23246, 12316, 23261, 5924,
                                                                       5942, 12696, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23741, 0, 3,
                                                                       23261, 12326, 23276, 5942,
                                                                       5960, 12726, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23786, 0, 3,
                                                                       23291, 12366, 23306, 5996,
                                                                       6014, 12816, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23831, 0, 3,
                                                                       23306, 12376, 23321, 6014,
                                                                       6032, 12846, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23876, 0, 3,
                                                                       23321, 12386, 23336, 6032,
                                                                       6050, 12876, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23921, 0, 3,
                                                                       23336, 12396, 23351, 6050,
                                                                       6068, 12906, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 23966, 0, 3,
                                                                       23351, 12406, 23366, 6068,
                                                                       6086, 12936, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24011, 0, 3,
                                                                       23366, 12416, 23381, 6086,
                                                                       6104, 12966, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24056, 0, 3,
                                                                       23381, 12426, 23396, 6104,
                                                                       6122, 12996, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 24101, 0, 3,
                                                                       23396, 12436, 23411, 6122,
                                                                       6140, 13026, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 24146, 0, 3,
                                                                       23426, 12516, 23471, 6176,
                                                                       6212, 13176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 24236, 0, 3,
                                                                       23471, 12546, 23516, 6212,
                                                                       6248, 13236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 24326, 0, 3,
                                                                       23516, 12576, 23561, 6248,
                                                                       6284, 13296, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 24416, 0, 3,
                                                                       23561, 12606, 23606, 6284,
                                                                       6320, 13356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 24506, 0, 3,
                                                                       23606, 12636, 23651, 6320,
                                                                       6356, 13416, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 24596, 0, 3,
                                                                       23651, 12666, 23696, 6356,
                                                                       6392, 13476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 24686, 0, 3,
                                                                       23696, 12696, 23741, 6392,
                                                                       6428, 13536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 24776, 0, 3,
                                                                       23786, 12816, 23831, 6500,
                                                                       6536, 13716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 24866, 0, 3,
                                                                       23831, 12846, 23876, 6536,
                                                                       6572, 13776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 24956, 0, 3,
                                                                       23876, 12876, 23921, 6572,
                                                                       6608, 13836, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25046, 0, 3,
                                                                       23921, 12906, 23966, 6608,
                                                                       6644, 13896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25136, 0, 3,
                                                                       23966, 12936, 24011, 6644,
                                                                       6680, 13956, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25226, 0, 3,
                                                                       24011, 12966, 24056, 6680,
                                                                       6716, 14016, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 25316, 0, 3,
                                                                       24056, 12996, 24101, 6716,
                                                                       6752, 14076, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 25406, 0, 3,
                                                                       24146, 13176, 24236, 6824,
                                                                       6884, 14336, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 25556, 0, 3,
                                                                       24236, 13236, 24326, 6884,
                                                                       6944, 14436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 25706, 0, 3,
                                                                       24326, 13296, 24416, 6944,
                                                                       7004, 14536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 25856, 0, 3,
                                                                       24416, 13356, 24506, 7004,
                                                                       7064, 14636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26006, 0, 3,
                                                                       24506, 13416, 24596, 7064,
                                                                       7124, 14736, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26156, 0, 3,
                                                                       24596, 13476, 24686, 7124,
                                                                       7184, 14836, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26306, 0, 3,
                                                                       24776, 13716, 24866, 7304,
                                                                       7364, 15136, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26456, 0, 3,
                                                                       24866, 13776, 24956, 7364,
                                                                       7424, 15236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26606, 0, 3,
                                                                       24956, 13836, 25046, 7424,
                                                                       7484, 15336, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26756, 0, 3,
                                                                       25046, 13896, 25136, 7484,
                                                                       7544, 15436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 26906, 0, 3,
                                                                       25136, 13956, 25226, 7544,
                                                                       7604, 15536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27056, 0, 3,
                                                                       25226, 14016, 25316, 7604,
                                                                       7664, 15636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27206, 0, 3,
                                                                       25406, 14336, 25556, 7784,
                                                                       7874, 16036, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27431, 0, 3,
                                                                       25556, 14436, 25706, 7874,
                                                                       7964, 16186, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27656, 0, 3,
                                                                       25706, 14536, 25856, 7964,
                                                                       8054, 16336, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 27881, 0, 3,
                                                                       25856, 14636, 26006, 8054,
                                                                       8144, 16486, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28106, 0, 3,
                                                                       26006, 14736, 26156, 8144,
                                                                       8234, 16636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28331, 0, 3,
                                                                       26306, 15136, 26456, 8414,
                                                                       8504, 17086, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28556, 0, 3,
                                                                       26456, 15236, 26606, 8504,
                                                                       8594, 17236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28781, 0, 3,
                                                                       26606, 15336, 26756, 8594,
                                                                       8684, 17386, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 29006, 0, 3,
                                                                       26756, 15436, 26906, 8684,
                                                                       8774, 17536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 29231, 0, 3,
                                                                       26906, 15536, 27056, 8774,
                                                                       8864, 17686, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29456, 0, 3,
                                                                       27206, 16036, 27431, 9044,
                                                                       9170, 18256, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 29771, 0, 3,
                                                                       27431, 16186, 27656, 9170,
                                                                       9296, 18466, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 30086, 0, 3,
                                                                       27656, 16336, 27881, 9296,
                                                                       9422, 18676, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 30401, 0, 3,
                                                                       27881, 16486, 28106, 9422,
                                                                       9548, 18886, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 30716, 0, 3,
                                                                       28331, 17086, 28556, 9800,
                                                                       9926, 19516, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 31031, 0, 3,
                                                                       28556, 17236, 28781, 9926,
                                                                       10052, 19726, ncols,
                                                                       gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 31346, 0, 3,
                                                                       28781, 17386, 29006,
                                                                       10052, 10178, 19936,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 31661, 0, 3,
                                                                       29006, 17536, 29231,
                                                                       10178, 10304, 20146,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 31976, 0, 3,
                                                                       29456, 18256, 29771,
                                                                       10556, 10724, 20916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 32396, 0, 3,
                                                                       29771, 18466, 30086,
                                                                       10724, 10892, 21196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 32816, 0, 3,
                                                                       30086, 18676, 30401,
                                                                       10892, 11060, 21476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 33236, 0, 3,
                                                                       30716, 19516, 31031,
                                                                       11396, 11564, 22316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 33656, 0, 3,
                                                                       31031, 19726, 31346,
                                                                       11564, 11732, 22596,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 34076, 0, 3,
                                                                       31346, 19936, 31661,
                                                                       11732, 11900, 22876,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34496, 3, 12236,
                                                                       12246, 23156, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34517, 3, 12246,
                                                                       12256, 23171, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34538, 3, 12256,
                                                                       12266, 23186, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34559, 3, 12266,
                                                                       12276, 23201, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34580, 3, 12276,
                                                                       12286, 23216, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34601, 3, 12286,
                                                                       12296, 23231, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34622, 3, 12296,
                                                                       12306, 23246, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34643, 3, 12306,
                                                                       12316, 23261, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34664, 3, 12316,
                                                                       12326, 23276, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34685, 3, 12346,
                                                                       12356, 23291, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34706, 3, 12356,
                                                                       12366, 23306, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34727, 3, 12366,
                                                                       12376, 23321, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34748, 3, 12376,
                                                                       12386, 23336, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34769, 3, 12386,
                                                                       12396, 23351, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34790, 3, 12396,
                                                                       12406, 23366, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34811, 3, 12406,
                                                                       12416, 23381, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34832, 3, 12416,
                                                                       12426, 23396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 34853, 3, 12426,
                                                                       12436, 23411, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 34874, 0, 3,
                                                                       34496, 23156, 34517,
                                                                       12456, 12486, 23426,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 34937, 0, 3,
                                                                       34517, 23171, 34538,
                                                                       12486, 12516, 23471,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35000, 0, 3,
                                                                       34538, 23186, 34559,
                                                                       12516, 12546, 23516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35063, 0, 3,
                                                                       34559, 23201, 34580,
                                                                       12546, 12576, 23561,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35126, 0, 3,
                                                                       34580, 23216, 34601,
                                                                       12576, 12606, 23606,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35189, 0, 3,
                                                                       34601, 23231, 34622,
                                                                       12606, 12636, 23651,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35252, 0, 3,
                                                                       34622, 23246, 34643,
                                                                       12636, 12666, 23696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35315, 0, 3,
                                                                       34643, 23261, 34664,
                                                                       12666, 12696, 23741,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35378, 0, 3,
                                                                       34685, 23291, 34706,
                                                                       12756, 12786, 23786,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35441, 0, 3,
                                                                       34706, 23306, 34727,
                                                                       12786, 12816, 23831,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35504, 0, 3,
                                                                       34727, 23321, 34748,
                                                                       12816, 12846, 23876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35567, 0, 3,
                                                                       34748, 23336, 34769,
                                                                       12846, 12876, 23921,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35630, 0, 3,
                                                                       34769, 23351, 34790,
                                                                       12876, 12906, 23966,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35693, 0, 3,
                                                                       34790, 23366, 34811,
                                                                       12906, 12936, 24011,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35756, 0, 3,
                                                                       34811, 23381, 34832,
                                                                       12936, 12966, 24056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35819, 0, 3,
                                                                       34832, 23396, 34853,
                                                                       12966, 12996, 24101,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 35882, 0, 3,
                                                                       34874, 23426, 34937,
                                                                       13056, 13116, 24146,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36008, 0, 3,
                                                                       34937, 23471, 35000,
                                                                       13116, 13176, 24236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36134, 0, 3,
                                                                       35000, 23516, 35063,
                                                                       13176, 13236, 24326,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36260, 0, 3,
                                                                       35063, 23561, 35126,
                                                                       13236, 13296, 24416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36386, 0, 3,
                                                                       35126, 23606, 35189,
                                                                       13296, 13356, 24506,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36512, 0, 3,
                                                                       35189, 23651, 35252,
                                                                       13356, 13416, 24596,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36638, 0, 3,
                                                                       35252, 23696, 35315,
                                                                       13416, 13476, 24686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36764, 0, 3,
                                                                       35378, 23786, 35441,
                                                                       13596, 13656, 24776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36890, 0, 3,
                                                                       35441, 23831, 35504,
                                                                       13656, 13716, 24866,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37016, 0, 3,
                                                                       35504, 23876, 35567,
                                                                       13716, 13776, 24956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37142, 0, 3,
                                                                       35567, 23921, 35630,
                                                                       13776, 13836, 25046,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37268, 0, 3,
                                                                       35630, 23966, 35693,
                                                                       13836, 13896, 25136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37394, 0, 3,
                                                                       35693, 24011, 35756,
                                                                       13896, 13956, 25226,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37520, 0, 3,
                                                                       35756, 24056, 35819,
                                                                       13956, 14016, 25316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37646, 0, 3,
                                                                       35882, 24146, 36008,
                                                                       14136, 14236, 25406,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37856, 0, 3,
                                                                       36008, 24236, 36134,
                                                                       14236, 14336, 25556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38066, 0, 3,
                                                                       36134, 24326, 36260,
                                                                       14336, 14436, 25706,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38276, 0, 3,
                                                                       36260, 24416, 36386,
                                                                       14436, 14536, 25856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38486, 0, 3,
                                                                       36386, 24506, 36512,
                                                                       14536, 14636, 26006,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38696, 0, 3,
                                                                       36512, 24596, 36638,
                                                                       14636, 14736, 26156,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38906, 0, 3,
                                                                       36764, 24776, 36890,
                                                                       14936, 15036, 26306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39116, 0, 3,
                                                                       36890, 24866, 37016,
                                                                       15036, 15136, 26456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39326, 0, 3,
                                                                       37016, 24956, 37142,
                                                                       15136, 15236, 26606,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39536, 0, 3,
                                                                       37142, 25046, 37268,
                                                                       15236, 15336, 26756,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39746, 0, 3,
                                                                       37268, 25136, 37394,
                                                                       15336, 15436, 26906,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 39956, 0, 3,
                                                                       37394, 25226, 37520,
                                                                       15436, 15536, 27056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40166, 0, 3,
                                                                       37646, 25406, 37856,
                                                                       15736, 15886, 27206,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40481, 0, 3,
                                                                       37856, 25556, 38066,
                                                                       15886, 16036, 27431,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 40796, 0, 3,
                                                                       38066, 25706, 38276,
                                                                       16036, 16186, 27656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41111, 0, 3,
                                                                       38276, 25856, 38486,
                                                                       16186, 16336, 27881,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41426, 0, 3,
                                                                       38486, 26006, 38696,
                                                                       16336, 16486, 28106,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 41741, 0, 3,
                                                                       38906, 26306, 39116,
                                                                       16786, 16936, 28331,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 42056, 0, 3,
                                                                       39116, 26456, 39326,
                                                                       16936, 17086, 28556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 42371, 0, 3,
                                                                       39326, 26606, 39536,
                                                                       17086, 17236, 28781,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 42686, 0, 3,
                                                                       39536, 26756, 39746,
                                                                       17236, 17386, 29006,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 43001, 0, 3,
                                                                       39746, 26906, 39956,
                                                                       17386, 17536, 29231,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 43316, 0, 3,
                                                                       40166, 27206, 40481,
                                                                       17836, 18046, 29456,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 43757, 0, 3,
                                                                       40481, 27431, 40796,
                                                                       18046, 18256, 29771,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 44198, 0, 3,
                                                                       40796, 27656, 41111,
                                                                       18256, 18466, 30086,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 44639, 0, 3,
                                                                       41111, 27881, 41426,
                                                                       18466, 18676, 30401,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 45080, 0, 3,
                                                                       41741, 28331, 42056,
                                                                       19096, 19306, 30716,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 45521, 0, 3,
                                                                       42056, 28556, 42371,
                                                                       19306, 19516, 31031,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 45962, 0, 3,
                                                                       42371, 28781, 42686,
                                                                       19516, 19726, 31346,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 46403, 0, 3,
                                                                       42686, 29006, 43001,
                                                                       19726, 19936, 31661,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 46844, 0, 3,
                                                                       43316, 29456, 43757,
                                                                       20356, 20636, 31976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 47432, 0, 3,
                                                                       43757, 29771, 44198,
                                                                       20636, 20916, 32396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 48020, 0, 3,
                                                                       44198, 30086, 44639,
                                                                       20916, 21196, 32816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 48608, 0, 3,
                                                                       45080, 30716, 45521,
                                                                       21756, 22036, 33236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 49196, 0, 3,
                                                                       45521, 31031, 45962,
                                                                       22036, 22316, 33656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 49784, 0, 3,
                                                                       45962, 31346, 46403,
                                                                       22316, 22596, 34076,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50372, 3, 23156,
                                                                       23171, 34538, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50400, 3, 23171,
                                                                       23186, 34559, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50428, 3, 23186,
                                                                       23201, 34580, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50456, 3, 23201,
                                                                       23216, 34601, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50484, 3, 23216,
                                                                       23231, 34622, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50512, 3, 23231,
                                                                       23246, 34643, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50540, 3, 23246,
                                                                       23261, 34664, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50568, 3, 23291,
                                                                       23306, 34727, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50596, 3, 23306,
                                                                       23321, 34748, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50624, 3, 23321,
                                                                       23336, 34769, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50652, 3, 23336,
                                                                       23351, 34790, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50680, 3, 23351,
                                                                       23366, 34811, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50708, 3, 23366,
                                                                       23381, 34832, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 50736, 3, 23381,
                                                                       23396, 34853, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 50764, 0, 3,
                                                                       50372, 34538, 50400,
                                                                       23426, 23471, 35000,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 50848, 0, 3,
                                                                       50400, 34559, 50428,
                                                                       23471, 23516, 35063,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 50932, 0, 3,
                                                                       50428, 34580, 50456,
                                                                       23516, 23561, 35126,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 51016, 0, 3,
                                                                       50456, 34601, 50484,
                                                                       23561, 23606, 35189,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 51100, 0, 3,
                                                                       50484, 34622, 50512,
                                                                       23606, 23651, 35252,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 51184, 0, 3,
                                                                       50512, 34643, 50540,
                                                                       23651, 23696, 35315,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 51268, 0, 3,
                                                                       50568, 34727, 50596,
                                                                       23786, 23831, 35504,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 51352, 0, 3,
                                                                       50596, 34748, 50624,
                                                                       23831, 23876, 35567,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 51436, 0, 3,
                                                                       50624, 34769, 50652,
                                                                       23876, 23921, 35630,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 51520, 0, 3,
                                                                       50652, 34790, 50680,
                                                                       23921, 23966, 35693,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 51604, 0, 3,
                                                                       50680, 34811, 50708,
                                                                       23966, 24011, 35756,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 51688, 0, 3,
                                                                       50708, 34832, 50736,
                                                                       24011, 24056, 35819,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 51772, 0, 3,
                                                                       50764, 35000, 50848,
                                                                       24146, 24236, 36134,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 51940, 0, 3,
                                                                       50848, 35063, 50932,
                                                                       24236, 24326, 36260,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 52108, 0, 3,
                                                                       50932, 35126, 51016,
                                                                       24326, 24416, 36386,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 52276, 0, 3,
                                                                       51016, 35189, 51100,
                                                                       24416, 24506, 36512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 52444, 0, 3,
                                                                       51100, 35252, 51184,
                                                                       24506, 24596, 36638,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 52612, 0, 3,
                                                                       51268, 35504, 51352,
                                                                       24776, 24866, 37016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 52780, 0, 3,
                                                                       51352, 35567, 51436,
                                                                       24866, 24956, 37142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 52948, 0, 3,
                                                                       51436, 35630, 51520,
                                                                       24956, 25046, 37268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 53116, 0, 3,
                                                                       51520, 35693, 51604,
                                                                       25046, 25136, 37394,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 53284, 0, 3,
                                                                       51604, 35756, 51688,
                                                                       25136, 25226, 37520,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 53452, 0, 3,
                                                                       51772, 36134, 51940,
                                                                       25406, 25556, 38066,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 53732, 0, 3,
                                                                       51940, 36260, 52108,
                                                                       25556, 25706, 38276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 54012, 0, 3,
                                                                       52108, 36386, 52276,
                                                                       25706, 25856, 38486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 54292, 0, 3,
                                                                       52276, 36512, 52444,
                                                                       25856, 26006, 38696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 54572, 0, 3,
                                                                       52612, 37016, 52780,
                                                                       26306, 26456, 39326,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 54852, 0, 3,
                                                                       52780, 37142, 52948,
                                                                       26456, 26606, 39536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 55132, 0, 3,
                                                                       52948, 37268, 53116,
                                                                       26606, 26756, 39746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 55412, 0, 3,
                                                                       53116, 37394, 53284,
                                                                       26756, 26906, 39956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 55692, 0, 3,
                                                                       53452, 38066, 53732,
                                                                       27206, 27431, 40796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 56112, 0, 3,
                                                                       53732, 38276, 54012,
                                                                       27431, 27656, 41111,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 56532, 0, 3,
                                                                       54012, 38486, 54292,
                                                                       27656, 27881, 41426,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 56952, 0, 3,
                                                                       54572, 39326, 54852,
                                                                       28331, 28556, 42371,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 57372, 0, 3,
                                                                       54852, 39536, 55132,
                                                                       28556, 28781, 42686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 57792, 0, 3,
                                                                       55132, 39746, 55412,
                                                                       28781, 29006, 43001,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 58212, 0, 3,
                                                                       55692, 40796, 56112,
                                                                       29456, 29771, 44198,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 58800, 0, 3,
                                                                       56112, 41111, 56532,
                                                                       29771, 30086, 44639,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 59388, 0, 3,
                                                                       56952, 42371, 57372,
                                                                       30716, 31031, 45962,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 59976, 0, 3,
                                                                       57372, 42686, 57792,
                                                                       31031, 31346, 46403,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 60564, 0, 3,
                                                                       58212, 44198, 58800,
                                                                       31976, 32396, 48020,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 61348, 0, 3,
                                                                       59388, 45962, 59976,
                                                                       33236, 33656, 49784,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62132, 3, 34496,
                                                                       34517, 50372, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62168, 3, 34517,
                                                                       34538, 50400, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62204, 3, 34538,
                                                                       34559, 50428, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62240, 3, 34559,
                                                                       34580, 50456, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62276, 3, 34580,
                                                                       34601, 50484, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62312, 3, 34601,
                                                                       34622, 50512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62348, 3, 34622,
                                                                       34643, 50540, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62384, 3, 34685,
                                                                       34706, 50568, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62420, 3, 34706,
                                                                       34727, 50596, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62456, 3, 34727,
                                                                       34748, 50624, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62492, 3, 34748,
                                                                       34769, 50652, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62528, 3, 34769,
                                                                       34790, 50680, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62564, 3, 34790,
                                                                       34811, 50708, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 62600, 3, 34811,
                                                                       34832, 50736, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 62636, 0, 3,
                                                                       62132, 50372, 62168,
                                                                       34874, 34937, 50764,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 62744, 0, 3,
                                                                       62168, 50400, 62204,
                                                                       34937, 35000, 50848,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 62852, 0, 3,
                                                                       62204, 50428, 62240,
                                                                       35000, 35063, 50932,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 62960, 0, 3,
                                                                       62240, 50456, 62276,
                                                                       35063, 35126, 51016,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 63068, 0, 3,
                                                                       62276, 50484, 62312,
                                                                       35126, 35189, 51100,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 63176, 0, 3,
                                                                       62312, 50512, 62348,
                                                                       35189, 35252, 51184,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 63284, 0, 3,
                                                                       62384, 50568, 62420,
                                                                       35378, 35441, 51268,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 63392, 0, 3,
                                                                       62420, 50596, 62456,
                                                                       35441, 35504, 51352,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 63500, 0, 3,
                                                                       62456, 50624, 62492,
                                                                       35504, 35567, 51436,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 63608, 0, 3,
                                                                       62492, 50652, 62528,
                                                                       35567, 35630, 51520,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 63716, 0, 3,
                                                                       62528, 50680, 62564,
                                                                       35630, 35693, 51604,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 63824, 0, 3,
                                                                       62564, 50708, 62600,
                                                                       35693, 35756, 51688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 63932, 0, 3,
                                                                       62636, 50764, 62744,
                                                                       35882, 36008, 51772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 64148, 0, 3,
                                                                       62744, 50848, 62852,
                                                                       36008, 36134, 51940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 64364, 0, 3,
                                                                       62852, 50932, 62960,
                                                                       36134, 36260, 52108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 64580, 0, 3,
                                                                       62960, 51016, 63068,
                                                                       36260, 36386, 52276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 64796, 0, 3,
                                                                       63068, 51100, 63176,
                                                                       36386, 36512, 52444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 65012, 0, 3,
                                                                       63284, 51268, 63392,
                                                                       36764, 36890, 52612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 65228, 0, 3,
                                                                       63392, 51352, 63500,
                                                                       36890, 37016, 52780,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 65444, 0, 3,
                                                                       63500, 51436, 63608,
                                                                       37016, 37142, 52948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 65660, 0, 3,
                                                                       63608, 51520, 63716,
                                                                       37142, 37268, 53116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 65876, 0, 3,
                                                                       63716, 51604, 63824,
                                                                       37268, 37394, 53284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 66092, 0, 3,
                                                                       63932, 51772, 64148,
                                                                       37646, 37856, 53452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 66452, 0, 3,
                                                                       64148, 51940, 64364,
                                                                       37856, 38066, 53732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 66812, 0, 3,
                                                                       64364, 52108, 64580,
                                                                       38066, 38276, 54012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 67172, 0, 3,
                                                                       64580, 52276, 64796,
                                                                       38276, 38486, 54292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 67532, 0, 3,
                                                                       65012, 52612, 65228,
                                                                       38906, 39116, 54572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 67892, 0, 3,
                                                                       65228, 52780, 65444,
                                                                       39116, 39326, 54852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 68252, 0, 3,
                                                                       65444, 52948, 65660,
                                                                       39326, 39536, 55132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 68612, 0, 3,
                                                                       65660, 53116, 65876,
                                                                       39536, 39746, 55412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 68972, 0, 3,
                                                                       66092, 53452, 66452,
                                                                       40166, 40481, 55692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 69512, 0, 3,
                                                                       66452, 53732, 66812,
                                                                       40481, 40796, 56112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 70052, 0, 3,
                                                                       66812, 54012, 67172,
                                                                       40796, 41111, 56532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 70592, 0, 3,
                                                                       67532, 54572, 67892,
                                                                       41741, 42056, 56952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 71132, 0, 3,
                                                                       67892, 54852, 68252,
                                                                       42056, 42371, 57372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 71672, 0, 3,
                                                                       68252, 55132, 68612,
                                                                       42371, 42686, 57792,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 72212, 0, 3,
                                                                       68972, 55692, 69512,
                                                                       43316, 43757, 58212,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 72968, 0, 3,
                                                                       69512, 56112, 70052,
                                                                       43757, 44198, 58800,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 73724, 0, 3,
                                                                       70592, 56952, 71132,
                                                                       45080, 45521, 59388,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 74480, 0, 3,
                                                                       71132, 57372, 71672,
                                                                       45521, 45962, 59976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 75236, 0, 3,
                                                                       72212, 58212, 72968,
                                                                       46844, 47432, 60564,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 76244, 0, 3,
                                                                       73724, 59388, 74480,
                                                                       48608, 49196, 61348,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 77252, 68972, 540, ncols);

                    simdfunc::contract_primitives(buffer, 78017, 70592, 540, ncols);

                    simdfunc::contract_primitives(buffer, 78782, 72212, 756, ncols);

                    simdfunc::contract_primitives(buffer, 79853, 73724, 756, ncols);

                    simdfunc::contract_primitives(buffer, 80924, 75236, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 82352, 76244, 1008, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 77792, 77252, 15, 1, nmax);

        simdtrf::transform_k_inner(buffer, 78557, 78017, 15, 1, nmax);

        simdtrf::transform_k_inner(buffer, 79538, 78782, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 80609, 79853, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 81932, 80924, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 83360, 82352, 28, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 83780, 77792, 79538, 15, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 84455, 78557, 80609, 15, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 85130, 79538, 81932, 15, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 86075, 80609, 83360, 15, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 87020, 83780, 85130, 15, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 88370, 84455, 86075, 15, nmax);

        simdtrf::transform_g_inner(buffer, 89720, 88370, 6, 15, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 89720, 135, nmax);

        simdtrf::transform_g_inner(buffer, 89720, 87020, 6, 15, nmax);

        simdtrf::transform_d_outer(values + 675 * nvalues + n * npairs, nvalues, buffer, 89720,
                                   135, nmax);
    }

    for (size_t m = 0; m < 1350; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
