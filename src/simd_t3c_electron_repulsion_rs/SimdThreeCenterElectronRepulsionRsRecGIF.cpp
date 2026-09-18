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


#include "SimdThreeCenterElectronRepulsionRsRecGIF.hpp"

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
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_gif_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_gif_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 84657, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1638 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 84657, 44588, 7358, dimensions);

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

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 708,
                                                                       729, 1044, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1472, 0, 3, 729,
                                                                       750, 1072, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1508, 0, 3, 750,
                                                                       771, 1100, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1544, 0, 3, 771,
                                                                       792, 1128, 1156, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1580, 0, 3, 792,
                                                                       813, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1616, 0, 3, 813,
                                                                       834, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1652, 0, 3, 876,
                                                                       897, 1240, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 897,
                                                                       918, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1724, 0, 3, 918,
                                                                       939, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1760, 0, 3, 939,
                                                                       960, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1796, 0, 3, 960,
                                                                       981, 1352, 1380, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1832, 0, 3, 981,
                                                                       1002, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1868, 0, 3, 1044,
                                                                       1072, 1436, 1472, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1913, 0, 3, 1072,
                                                                       1100, 1472, 1508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1958, 0, 3, 1100,
                                                                       1128, 1508, 1544, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2003, 0, 3, 1128,
                                                                       1156, 1544, 1580, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2048, 0, 3, 1156,
                                                                       1184, 1580, 1616, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2093, 0, 3, 1240,
                                                                       1268, 1652, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2138, 0, 3, 1268,
                                                                       1296, 1688, 1724, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2183, 0, 3, 1296,
                                                                       1324, 1724, 1760, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2228, 0, 3, 1324,
                                                                       1352, 1760, 1796, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2273, 0, 3, 1352,
                                                                       1380, 1796, 1832, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2318, 0, 3, 1436,
                                                                       1472, 1868, 1913, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2373, 0, 3, 1472,
                                                                       1508, 1913, 1958, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2428, 0, 3, 1508,
                                                                       1544, 1958, 2003, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2483, 0, 3, 1544,
                                                                       1580, 2003, 2048, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2538, 0, 3, 1652,
                                                                       1688, 2093, 2138, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2593, 0, 3, 1688,
                                                                       1724, 2138, 2183, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1724,
                                                                       1760, 2183, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2703, 0, 3, 1760,
                                                                       1796, 2228, 2273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2758, 0, 3, 1868,
                                                                       1913, 2318, 2373, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2824, 0, 3, 1913,
                                                                       1958, 2373, 2428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2890, 0, 3, 1958,
                                                                       2003, 2428, 2483, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2956, 0, 3, 2093,
                                                                       2138, 2538, 2593, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3022, 0, 3, 2138,
                                                                       2183, 2593, 2648, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 2183,
                                                                       2228, 2648, 2703, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3154, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3157, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3160, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3163, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3166, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3169, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3172, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3175, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3178, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3181, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3184, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3187, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3190, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3193, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3196, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3199, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3202, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3205, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3208, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3211, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3214, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3217, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3220, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3223, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3226, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3229, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3232, 3, 9, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3241, 3, 10, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3250, 3, 11, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3259, 3, 12, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3268, 3, 13, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3277, 3, 14, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3286, 3, 15, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3295, 3, 16, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3304, 3, 17, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3313, 3, 18, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3322, 3, 23, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3331, 3, 24, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3340, 3, 25, 82,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3349, 3, 26, 85,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3358, 3, 27, 88,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3367, 3, 28, 91,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3376, 3, 29, 94,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3385, 3, 30, 97,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3394, 3, 31, 100,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3403, 3, 32, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3412, 3, 34, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3430, 3, 37, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3448, 3, 40, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3466, 3, 43, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3484, 3, 46, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3502, 3, 49, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3520, 3, 52, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3538, 3, 55, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3556, 3, 58, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3574, 3, 61, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3592, 3, 64, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3610, 3, 70, 172,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3628, 3, 73, 178,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3646, 3, 76, 184,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3664, 3, 79, 190,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3682, 3, 82, 196,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3700, 3, 85, 202,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3718, 3, 88, 208,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3736, 3, 91, 214,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3754, 3, 94, 220,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3772, 3, 97, 226,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3790, 3, 100, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3808, 3, 106, 238,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3838, 3, 112, 248,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3868, 3, 118, 258,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3898, 3, 124, 268,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3928, 3, 130, 278,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3958, 3, 136, 288,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3988, 3, 142, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4018, 3, 148, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4048, 3, 154, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4078, 3, 160, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4108, 3, 172, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4138, 3, 178, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4168, 3, 184, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4198, 3, 190, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4228, 3, 196, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4258, 3, 202, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4288, 3, 208, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4318, 3, 214, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4348, 3, 220, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4378, 3, 226, 428,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4408, 3, 238, 438,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4453, 3, 248, 453,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4498, 3, 258, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4543, 3, 268, 483,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4588, 3, 278, 498,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4633, 3, 288, 513,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4678, 3, 298, 528,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4723, 3, 308, 543,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4768, 3, 318, 558,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4813, 3, 338, 573,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4858, 3, 348, 588,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4903, 3, 358, 603,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4948, 3, 368, 618,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4993, 3, 378, 633,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5038, 3, 388, 648,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5083, 3, 398, 663,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5128, 3, 408, 678,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5173, 3, 418, 693,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5218, 3, 438, 708,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5281, 3, 453, 729,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5344, 3, 468, 750,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5407, 3, 483, 771,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5470, 3, 498, 792,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5533, 3, 513, 813,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5596, 3, 528, 834,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5659, 3, 543, 855,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5722, 3, 573, 876,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5785, 3, 588, 897,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5848, 3, 603, 918,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5911, 3, 618, 939,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5974, 3, 633, 960,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6037, 3, 648, 981,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6100, 3, 663,
                                                                       1002, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6163, 3, 678,
                                                                       1023, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6226, 3, 708,
                                                                       1044, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6310, 3, 729,
                                                                       1072, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6394, 3, 750,
                                                                       1100, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6478, 3, 771,
                                                                       1128, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6562, 3, 792,
                                                                       1156, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6646, 3, 813,
                                                                       1184, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6730, 3, 834,
                                                                       1212, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6814, 3, 876,
                                                                       1240, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6898, 3, 897,
                                                                       1268, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6982, 3, 918,
                                                                       1296, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7066, 3, 939,
                                                                       1324, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7150, 3, 960,
                                                                       1352, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7234, 3, 981,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7318, 3, 1002,
                                                                       1408, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7402, 3, 1044,
                                                                       1436, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7510, 3, 1072,
                                                                       1472, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7618, 3, 1100,
                                                                       1508, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7726, 3, 1128,
                                                                       1544, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7834, 3, 1156,
                                                                       1580, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7942, 3, 1184,
                                                                       1616, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8050, 3, 1240,
                                                                       1652, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8158, 3, 1268,
                                                                       1688, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8266, 3, 1296,
                                                                       1724, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8374, 3, 1324,
                                                                       1760, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8482, 3, 1352,
                                                                       1796, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8590, 3, 1380,
                                                                       1832, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8698, 3, 1436,
                                                                       1868, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8833, 3, 1472,
                                                                       1913, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8968, 3, 1508,
                                                                       1958, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9103, 3, 1544,
                                                                       2003, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9238, 3, 1580,
                                                                       2048, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9373, 3, 1652,
                                                                       2093, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9508, 3, 1688,
                                                                       2138, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9643, 3, 1724,
                                                                       2183, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9778, 3, 1760,
                                                                       2228, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9913, 3, 1796,
                                                                       2273, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10048, 3, 1868,
                                                                       2318, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10213, 3, 1913,
                                                                       2373, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10378, 3, 1958,
                                                                       2428, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10543, 3, 2003,
                                                                       2483, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10708, 3, 2093,
                                                                       2538, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10873, 3, 2138,
                                                                       2593, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11038, 3, 2183,
                                                                       2648, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11203, 3, 2228,
                                                                       2703, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11368, 3, 2318,
                                                                       2758, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11566, 3, 2373,
                                                                       2824, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11764, 3, 2428,
                                                                       2890, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11962, 3, 2538,
                                                                       2956, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12160, 3, 2593,
                                                                       3022, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12358, 3, 2648,
                                                                       3088, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12556, 3, 7, 8,
                                                                       3160, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12562, 3, 8, 9,
                                                                       3163, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12568, 3, 9, 10,
                                                                       3166, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12574, 3, 10, 11,
                                                                       3169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12580, 3, 11, 12,
                                                                       3172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12586, 3, 12, 13,
                                                                       3175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12592, 3, 13, 14,
                                                                       3178, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12598, 3, 14, 15,
                                                                       3181, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12604, 3, 15, 16,
                                                                       3184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12610, 3, 16, 17,
                                                                       3187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12616, 3, 17, 18,
                                                                       3190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12622, 3, 21, 22,
                                                                       3199, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12628, 3, 22, 23,
                                                                       3202, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12634, 3, 23, 24,
                                                                       3205, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12640, 3, 24, 25,
                                                                       3208, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12646, 3, 25, 26,
                                                                       3211, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12652, 3, 26, 27,
                                                                       3214, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12658, 3, 27, 28,
                                                                       3217, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12664, 3, 28, 29,
                                                                       3220, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12670, 3, 29, 30,
                                                                       3223, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12676, 3, 30, 31,
                                                                       3226, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12682, 3, 31, 32,
                                                                       3229, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12688, 0, 3,
                                                                       12556, 3160, 12562, 3232,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12706, 0, 3,
                                                                       12562, 3163, 12568, 3241,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12724, 0, 3,
                                                                       12568, 3166, 12574, 3250,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12742, 0, 3,
                                                                       12574, 3169, 12580, 3259,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12760, 0, 3,
                                                                       12580, 3172, 12586, 3268,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12778, 0, 3,
                                                                       12586, 3175, 12592, 3277,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12796, 0, 3,
                                                                       12592, 3178, 12598, 3286,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12814, 0, 3,
                                                                       12598, 3181, 12604, 3295,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12832, 0, 3,
                                                                       12604, 3184, 12610, 3304,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12850, 0, 3,
                                                                       12610, 3187, 12616, 3313,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12868, 0, 3,
                                                                       12622, 3199, 12628, 3322,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12886, 0, 3,
                                                                       12628, 3202, 12634, 3331,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12904, 0, 3,
                                                                       12634, 3205, 12640, 3340,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12922, 0, 3,
                                                                       12640, 3208, 12646, 3349,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12940, 0, 3,
                                                                       12646, 3211, 12652, 3358,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12958, 0, 3,
                                                                       12652, 3214, 12658, 3367,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12976, 0, 3,
                                                                       12658, 3217, 12664, 3376,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12994, 0, 3,
                                                                       12664, 3220, 12670, 3385,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13012, 0, 3,
                                                                       12670, 3223, 12676, 3394,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 13030, 0, 3,
                                                                       12676, 3226, 12682, 3403,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13048, 0, 3,
                                                                       12688, 3232, 12706, 106,
                                                                       112, 3448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13084, 0, 3,
                                                                       12706, 3241, 12724, 112,
                                                                       118, 3466, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13120, 0, 3,
                                                                       12724, 3250, 12742, 118,
                                                                       124, 3484, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13156, 0, 3,
                                                                       12742, 3259, 12760, 124,
                                                                       130, 3502, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13192, 0, 3,
                                                                       12760, 3268, 12778, 130,
                                                                       136, 3520, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13228, 0, 3,
                                                                       12778, 3277, 12796, 136,
                                                                       142, 3538, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13264, 0, 3,
                                                                       12796, 3286, 12814, 142,
                                                                       148, 3556, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13300, 0, 3,
                                                                       12814, 3295, 12832, 148,
                                                                       154, 3574, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13336, 0, 3,
                                                                       12832, 3304, 12850, 154,
                                                                       160, 3592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13372, 0, 3,
                                                                       12868, 3322, 12886, 172,
                                                                       178, 3646, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13408, 0, 3,
                                                                       12886, 3331, 12904, 178,
                                                                       184, 3664, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13444, 0, 3,
                                                                       12904, 3340, 12922, 184,
                                                                       190, 3682, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13480, 0, 3,
                                                                       12922, 3349, 12940, 190,
                                                                       196, 3700, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13516, 0, 3,
                                                                       12940, 3358, 12958, 196,
                                                                       202, 3718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13552, 0, 3,
                                                                       12958, 3367, 12976, 202,
                                                                       208, 3736, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13588, 0, 3,
                                                                       12976, 3376, 12994, 208,
                                                                       214, 3754, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13624, 0, 3,
                                                                       12994, 3385, 13012, 214,
                                                                       220, 3772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 13660, 0, 3,
                                                                       13012, 3394, 13030, 220,
                                                                       226, 3790, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13696, 0, 3,
                                                                       13048, 3448, 13084, 238,
                                                                       248, 3868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13756, 0, 3,
                                                                       13084, 3466, 13120, 248,
                                                                       258, 3898, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13816, 0, 3,
                                                                       13120, 3484, 13156, 258,
                                                                       268, 3928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13876, 0, 3,
                                                                       13156, 3502, 13192, 268,
                                                                       278, 3958, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13936, 0, 3,
                                                                       13192, 3520, 13228, 278,
                                                                       288, 3988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13996, 0, 3,
                                                                       13228, 3538, 13264, 288,
                                                                       298, 4018, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14056, 0, 3,
                                                                       13264, 3556, 13300, 298,
                                                                       308, 4048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14116, 0, 3,
                                                                       13300, 3574, 13336, 308,
                                                                       318, 4078, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14176, 0, 3,
                                                                       13372, 3646, 13408, 338,
                                                                       348, 4168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14236, 0, 3,
                                                                       13408, 3664, 13444, 348,
                                                                       358, 4198, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14296, 0, 3,
                                                                       13444, 3682, 13480, 358,
                                                                       368, 4228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14356, 0, 3,
                                                                       13480, 3700, 13516, 368,
                                                                       378, 4258, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14416, 0, 3,
                                                                       13516, 3718, 13552, 378,
                                                                       388, 4288, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14476, 0, 3,
                                                                       13552, 3736, 13588, 388,
                                                                       398, 4318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14536, 0, 3,
                                                                       13588, 3754, 13624, 398,
                                                                       408, 4348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 14596, 0, 3,
                                                                       13624, 3772, 13660, 408,
                                                                       418, 4378, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14656, 0, 3,
                                                                       13696, 3868, 13756, 438,
                                                                       453, 4498, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14746, 0, 3,
                                                                       13756, 3898, 13816, 453,
                                                                       468, 4543, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14836, 0, 3,
                                                                       13816, 3928, 13876, 468,
                                                                       483, 4588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14926, 0, 3,
                                                                       13876, 3958, 13936, 483,
                                                                       498, 4633, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15016, 0, 3,
                                                                       13936, 3988, 13996, 498,
                                                                       513, 4678, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15106, 0, 3,
                                                                       13996, 4018, 14056, 513,
                                                                       528, 4723, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15196, 0, 3,
                                                                       14056, 4048, 14116, 528,
                                                                       543, 4768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15286, 0, 3,
                                                                       14176, 4168, 14236, 573,
                                                                       588, 4903, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15376, 0, 3,
                                                                       14236, 4198, 14296, 588,
                                                                       603, 4948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15466, 0, 3,
                                                                       14296, 4228, 14356, 603,
                                                                       618, 4993, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15556, 0, 3,
                                                                       14356, 4258, 14416, 618,
                                                                       633, 5038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15646, 0, 3,
                                                                       14416, 4288, 14476, 633,
                                                                       648, 5083, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15736, 0, 3,
                                                                       14476, 4318, 14536, 648,
                                                                       663, 5128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15826, 0, 3,
                                                                       14536, 4348, 14596, 663,
                                                                       678, 5173, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15916, 0, 3,
                                                                       14656, 4498, 14746, 708,
                                                                       729, 5344, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16042, 0, 3,
                                                                       14746, 4543, 14836, 729,
                                                                       750, 5407, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16168, 0, 3,
                                                                       14836, 4588, 14926, 750,
                                                                       771, 5470, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16294, 0, 3,
                                                                       14926, 4633, 15016, 771,
                                                                       792, 5533, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16420, 0, 3,
                                                                       15016, 4678, 15106, 792,
                                                                       813, 5596, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16546, 0, 3,
                                                                       15106, 4723, 15196, 813,
                                                                       834, 5659, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16672, 0, 3,
                                                                       15286, 4903, 15376, 876,
                                                                       897, 5848, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16798, 0, 3,
                                                                       15376, 4948, 15466, 897,
                                                                       918, 5911, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16924, 0, 3,
                                                                       15466, 4993, 15556, 918,
                                                                       939, 5974, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17050, 0, 3,
                                                                       15556, 5038, 15646, 939,
                                                                       960, 6037, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17176, 0, 3,
                                                                       15646, 5083, 15736, 960,
                                                                       981, 6100, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17302, 0, 3,
                                                                       15736, 5128, 15826, 981,
                                                                       1002, 6163, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17428, 0, 3,
                                                                       15916, 5344, 16042, 1044,
                                                                       1072, 6394, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17596, 0, 3,
                                                                       16042, 5407, 16168, 1072,
                                                                       1100, 6478, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17764, 0, 3,
                                                                       16168, 5470, 16294, 1100,
                                                                       1128, 6562, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17932, 0, 3,
                                                                       16294, 5533, 16420, 1128,
                                                                       1156, 6646, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18100, 0, 3,
                                                                       16420, 5596, 16546, 1156,
                                                                       1184, 6730, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18268, 0, 3,
                                                                       16672, 5848, 16798, 1240,
                                                                       1268, 6982, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18436, 0, 3,
                                                                       16798, 5911, 16924, 1268,
                                                                       1296, 7066, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18604, 0, 3,
                                                                       16924, 5974, 17050, 1296,
                                                                       1324, 7150, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18772, 0, 3,
                                                                       17050, 6037, 17176, 1324,
                                                                       1352, 7234, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18940, 0, 3,
                                                                       17176, 6100, 17302, 1352,
                                                                       1380, 7318, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19108, 0, 3,
                                                                       17428, 6394, 17596, 1436,
                                                                       1472, 7618, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19324, 0, 3,
                                                                       17596, 6478, 17764, 1472,
                                                                       1508, 7726, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19540, 0, 3,
                                                                       17764, 6562, 17932, 1508,
                                                                       1544, 7834, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19756, 0, 3,
                                                                       17932, 6646, 18100, 1544,
                                                                       1580, 7942, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19972, 0, 3,
                                                                       18268, 6982, 18436, 1652,
                                                                       1688, 8266, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20188, 0, 3,
                                                                       18436, 7066, 18604, 1688,
                                                                       1724, 8374, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20404, 0, 3,
                                                                       18604, 7150, 18772, 1724,
                                                                       1760, 8482, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20620, 0, 3,
                                                                       18772, 7234, 18940, 1760,
                                                                       1796, 8590, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 20836, 0, 3,
                                                                       19108, 7618, 19324, 1868,
                                                                       1913, 8968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21106, 0, 3,
                                                                       19324, 7726, 19540, 1913,
                                                                       1958, 9103, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21376, 0, 3,
                                                                       19540, 7834, 19756, 1958,
                                                                       2003, 9238, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21646, 0, 3,
                                                                       19972, 8266, 20188, 2093,
                                                                       2138, 9643, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21916, 0, 3,
                                                                       20188, 8374, 20404, 2138,
                                                                       2183, 9778, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 22186, 0, 3,
                                                                       20404, 8482, 20620, 2183,
                                                                       2228, 9913, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 22456, 0, 3,
                                                                       20836, 8968, 21106, 2318,
                                                                       2373, 10378, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 22786, 0, 3,
                                                                       21106, 9103, 21376, 2373,
                                                                       2428, 10543, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 23116, 0, 3,
                                                                       21646, 9643, 21916, 2538,
                                                                       2593, 11038, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 23446, 0, 3,
                                                                       21916, 9778, 22186, 2593,
                                                                       2648, 11203, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 23776, 0, 3,
                                                                       22456, 10378, 22786, 2758,
                                                                       2824, 11764, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 24172, 0, 3,
                                                                       23116, 11038, 23446, 2956,
                                                                       3022, 12358, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24568, 3, 3154,
                                                                       3157, 12556, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24578, 3, 3157,
                                                                       3160, 12562, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24588, 3, 3160,
                                                                       3163, 12568, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24598, 3, 3163,
                                                                       3166, 12574, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24608, 3, 3166,
                                                                       3169, 12580, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24618, 3, 3169,
                                                                       3172, 12586, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24628, 3, 3172,
                                                                       3175, 12592, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24638, 3, 3175,
                                                                       3178, 12598, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24648, 3, 3178,
                                                                       3181, 12604, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24658, 3, 3181,
                                                                       3184, 12610, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24668, 3, 3184,
                                                                       3187, 12616, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24678, 3, 3193,
                                                                       3196, 12622, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24688, 3, 3196,
                                                                       3199, 12628, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24698, 3, 3199,
                                                                       3202, 12634, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24708, 3, 3202,
                                                                       3205, 12640, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24718, 3, 3205,
                                                                       3208, 12646, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24728, 3, 3208,
                                                                       3211, 12652, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24738, 3, 3211,
                                                                       3214, 12658, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24748, 3, 3214,
                                                                       3217, 12664, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24758, 3, 3217,
                                                                       3220, 12670, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24768, 3, 3220,
                                                                       3223, 12676, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 24778, 3, 3223,
                                                                       3226, 12682, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24788, 0, 3,
                                                                       24568, 12556, 24578,
                                                                       12688, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24818, 0, 3,
                                                                       24578, 12562, 24588,
                                                                       12706, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24848, 0, 3,
                                                                       24588, 12568, 24598,
                                                                       12724, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24878, 0, 3,
                                                                       24598, 12574, 24608,
                                                                       12742, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24908, 0, 3,
                                                                       24608, 12580, 24618,
                                                                       12760, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24938, 0, 3,
                                                                       24618, 12586, 24628,
                                                                       12778, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24968, 0, 3,
                                                                       24628, 12592, 24638,
                                                                       12796, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24998, 0, 3,
                                                                       24638, 12598, 24648,
                                                                       12814, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25028, 0, 3,
                                                                       24648, 12604, 24658,
                                                                       12832, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25058, 0, 3,
                                                                       24658, 12610, 24668,
                                                                       12850, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25088, 0, 3,
                                                                       24678, 12622, 24688,
                                                                       12868, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25118, 0, 3,
                                                                       24688, 12628, 24698,
                                                                       12886, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25148, 0, 3,
                                                                       24698, 12634, 24708,
                                                                       12904, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25178, 0, 3,
                                                                       24708, 12640, 24718,
                                                                       12922, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25208, 0, 3,
                                                                       24718, 12646, 24728,
                                                                       12940, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25238, 0, 3,
                                                                       24728, 12652, 24738,
                                                                       12958, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25268, 0, 3,
                                                                       24738, 12658, 24748,
                                                                       12976, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25298, 0, 3,
                                                                       24748, 12664, 24758,
                                                                       12994, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25328, 0, 3,
                                                                       24758, 12670, 24768,
                                                                       13012, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 25358, 0, 3,
                                                                       24768, 12676, 24778,
                                                                       13030, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25388, 0, 3,
                                                                       24788, 12688, 24818, 3412,
                                                                       3430, 13048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25448, 0, 3,
                                                                       24818, 12706, 24848, 3430,
                                                                       3448, 13084, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25508, 0, 3,
                                                                       24848, 12724, 24878, 3448,
                                                                       3466, 13120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25568, 0, 3,
                                                                       24878, 12742, 24908, 3466,
                                                                       3484, 13156, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25628, 0, 3,
                                                                       24908, 12760, 24938, 3484,
                                                                       3502, 13192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25688, 0, 3,
                                                                       24938, 12778, 24968, 3502,
                                                                       3520, 13228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25748, 0, 3,
                                                                       24968, 12796, 24998, 3520,
                                                                       3538, 13264, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25808, 0, 3,
                                                                       24998, 12814, 25028, 3538,
                                                                       3556, 13300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25868, 0, 3,
                                                                       25028, 12832, 25058, 3556,
                                                                       3574, 13336, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25928, 0, 3,
                                                                       25088, 12868, 25118, 3610,
                                                                       3628, 13372, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25988, 0, 3,
                                                                       25118, 12886, 25148, 3628,
                                                                       3646, 13408, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26048, 0, 3,
                                                                       25148, 12904, 25178, 3646,
                                                                       3664, 13444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26108, 0, 3,
                                                                       25178, 12922, 25208, 3664,
                                                                       3682, 13480, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26168, 0, 3,
                                                                       25208, 12940, 25238, 3682,
                                                                       3700, 13516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26228, 0, 3,
                                                                       25238, 12958, 25268, 3700,
                                                                       3718, 13552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26288, 0, 3,
                                                                       25268, 12976, 25298, 3718,
                                                                       3736, 13588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26348, 0, 3,
                                                                       25298, 12994, 25328, 3736,
                                                                       3754, 13624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 26408, 0, 3,
                                                                       25328, 13012, 25358, 3754,
                                                                       3772, 13660, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26468, 0, 3,
                                                                       25388, 13048, 25448, 3808,
                                                                       3838, 13696, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26568, 0, 3,
                                                                       25448, 13084, 25508, 3838,
                                                                       3868, 13756, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26668, 0, 3,
                                                                       25508, 13120, 25568, 3868,
                                                                       3898, 13816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26768, 0, 3,
                                                                       25568, 13156, 25628, 3898,
                                                                       3928, 13876, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26868, 0, 3,
                                                                       25628, 13192, 25688, 3928,
                                                                       3958, 13936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26968, 0, 3,
                                                                       25688, 13228, 25748, 3958,
                                                                       3988, 13996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 27068, 0, 3,
                                                                       25748, 13264, 25808, 3988,
                                                                       4018, 14056, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 27168, 0, 3,
                                                                       25808, 13300, 25868, 4018,
                                                                       4048, 14116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 27268, 0, 3,
                                                                       25928, 13372, 25988, 4108,
                                                                       4138, 14176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 27368, 0, 3,
                                                                       25988, 13408, 26048, 4138,
                                                                       4168, 14236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 27468, 0, 3,
                                                                       26048, 13444, 26108, 4168,
                                                                       4198, 14296, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 27568, 0, 3,
                                                                       26108, 13480, 26168, 4198,
                                                                       4228, 14356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 27668, 0, 3,
                                                                       26168, 13516, 26228, 4228,
                                                                       4258, 14416, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 27768, 0, 3,
                                                                       26228, 13552, 26288, 4258,
                                                                       4288, 14476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 27868, 0, 3,
                                                                       26288, 13588, 26348, 4288,
                                                                       4318, 14536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 27968, 0, 3,
                                                                       26348, 13624, 26408, 4318,
                                                                       4348, 14596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28068, 0, 3,
                                                                       26468, 13696, 26568, 4408,
                                                                       4453, 14656, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28218, 0, 3,
                                                                       26568, 13756, 26668, 4453,
                                                                       4498, 14746, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28368, 0, 3,
                                                                       26668, 13816, 26768, 4498,
                                                                       4543, 14836, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28518, 0, 3,
                                                                       26768, 13876, 26868, 4543,
                                                                       4588, 14926, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28668, 0, 3,
                                                                       26868, 13936, 26968, 4588,
                                                                       4633, 15016, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28818, 0, 3,
                                                                       26968, 13996, 27068, 4633,
                                                                       4678, 15106, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 28968, 0, 3,
                                                                       27068, 14056, 27168, 4678,
                                                                       4723, 15196, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 29118, 0, 3,
                                                                       27268, 14176, 27368, 4813,
                                                                       4858, 15286, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 29268, 0, 3,
                                                                       27368, 14236, 27468, 4858,
                                                                       4903, 15376, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 29418, 0, 3,
                                                                       27468, 14296, 27568, 4903,
                                                                       4948, 15466, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 29568, 0, 3,
                                                                       27568, 14356, 27668, 4948,
                                                                       4993, 15556, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 29718, 0, 3,
                                                                       27668, 14416, 27768, 4993,
                                                                       5038, 15646, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 29868, 0, 3,
                                                                       27768, 14476, 27868, 5038,
                                                                       5083, 15736, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 30018, 0, 3,
                                                                       27868, 14536, 27968, 5083,
                                                                       5128, 15826, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30168, 0, 3,
                                                                       28068, 14656, 28218, 5218,
                                                                       5281, 15916, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30378, 0, 3,
                                                                       28218, 14746, 28368, 5281,
                                                                       5344, 16042, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30588, 0, 3,
                                                                       28368, 14836, 28518, 5344,
                                                                       5407, 16168, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 30798, 0, 3,
                                                                       28518, 14926, 28668, 5407,
                                                                       5470, 16294, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 31008, 0, 3,
                                                                       28668, 15016, 28818, 5470,
                                                                       5533, 16420, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 31218, 0, 3,
                                                                       28818, 15106, 28968, 5533,
                                                                       5596, 16546, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 31428, 0, 3,
                                                                       29118, 15286, 29268, 5722,
                                                                       5785, 16672, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 31638, 0, 3,
                                                                       29268, 15376, 29418, 5785,
                                                                       5848, 16798, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 31848, 0, 3,
                                                                       29418, 15466, 29568, 5848,
                                                                       5911, 16924, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 32058, 0, 3,
                                                                       29568, 15556, 29718, 5911,
                                                                       5974, 17050, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 32268, 0, 3,
                                                                       29718, 15646, 29868, 5974,
                                                                       6037, 17176, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 32478, 0, 3,
                                                                       29868, 15736, 30018, 6037,
                                                                       6100, 17302, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32688, 0, 3,
                                                                       30168, 15916, 30378, 6226,
                                                                       6310, 17428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32968, 0, 3,
                                                                       30378, 16042, 30588, 6310,
                                                                       6394, 17596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 33248, 0, 3,
                                                                       30588, 16168, 30798, 6394,
                                                                       6478, 17764, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 33528, 0, 3,
                                                                       30798, 16294, 31008, 6478,
                                                                       6562, 17932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 33808, 0, 3,
                                                                       31008, 16420, 31218, 6562,
                                                                       6646, 18100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 34088, 0, 3,
                                                                       31428, 16672, 31638, 6814,
                                                                       6898, 18268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 34368, 0, 3,
                                                                       31638, 16798, 31848, 6898,
                                                                       6982, 18436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 34648, 0, 3,
                                                                       31848, 16924, 32058, 6982,
                                                                       7066, 18604, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 34928, 0, 3,
                                                                       32058, 17050, 32268, 7066,
                                                                       7150, 18772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 35208, 0, 3,
                                                                       32268, 17176, 32478, 7150,
                                                                       7234, 18940, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 35488, 0, 3,
                                                                       32688, 17428, 32968, 7402,
                                                                       7510, 19108, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 35848, 0, 3,
                                                                       32968, 17596, 33248, 7510,
                                                                       7618, 19324, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 36208, 0, 3,
                                                                       33248, 17764, 33528, 7618,
                                                                       7726, 19540, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 36568, 0, 3,
                                                                       33528, 17932, 33808, 7726,
                                                                       7834, 19756, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 36928, 0, 3,
                                                                       34088, 18268, 34368, 8050,
                                                                       8158, 19972, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 37288, 0, 3,
                                                                       34368, 18436, 34648, 8158,
                                                                       8266, 20188, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 37648, 0, 3,
                                                                       34648, 18604, 34928, 8266,
                                                                       8374, 20404, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 38008, 0, 3,
                                                                       34928, 18772, 35208, 8374,
                                                                       8482, 20620, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 38368, 0, 3,
                                                                       35488, 19108, 35848, 8698,
                                                                       8833, 20836, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 38818, 0, 3,
                                                                       35848, 19324, 36208, 8833,
                                                                       8968, 21106, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 39268, 0, 3,
                                                                       36208, 19540, 36568, 8968,
                                                                       9103, 21376, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 39718, 0, 3,
                                                                       36928, 19972, 37288, 9373,
                                                                       9508, 21646, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 40168, 0, 3,
                                                                       37288, 20188, 37648, 9508,
                                                                       9643, 21916, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 40618, 0, 3,
                                                                       37648, 20404, 38008, 9643,
                                                                       9778, 22186, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 41068, 0, 3,
                                                                       38368, 20836, 38818,
                                                                       10048, 10213, 22456,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 41618, 0, 3,
                                                                       38818, 21106, 39268,
                                                                       10213, 10378, 22786,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 42168, 0, 3,
                                                                       39718, 21646, 40168,
                                                                       10708, 10873, 23116,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 42718, 0, 3,
                                                                       40168, 21916, 40618,
                                                                       10873, 11038, 23446,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 43268, 0, 3,
                                                                       41068, 22456, 41618,
                                                                       11368, 11566, 23776,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 43928, 0, 3,
                                                                       42168, 23116, 42718,
                                                                       11962, 12160, 24172,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 44588, 32688, 280, ncols);

                    simdfunc::contract_primitives(buffer, 45064, 34088, 280, ncols);

                    simdfunc::contract_primitives(buffer, 45540, 35488, 360, ncols);

                    simdfunc::contract_primitives(buffer, 46152, 36928, 360, ncols);

                    simdfunc::contract_primitives(buffer, 46764, 38368, 450, ncols);

                    simdfunc::contract_primitives(buffer, 47529, 39718, 450, ncols);

                    simdfunc::contract_primitives(buffer, 48294, 41068, 550, ncols);

                    simdfunc::contract_primitives(buffer, 49229, 42168, 550, ncols);

                    simdfunc::contract_primitives(buffer, 50164, 43268, 660, ncols);

                    simdfunc::contract_primitives(buffer, 51286, 43928, 660, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 44868, 44588, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 45344, 45064, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 45900, 45540, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 46512, 46152, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 47214, 46764, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 47979, 47529, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 48844, 48294, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 49779, 49229, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 50824, 50164, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 51946, 51286, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 52408, 44868, 45900, 7, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 52996, 45344, 46512, 7, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 53584, 45900, 47214, 7, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 54340, 46512, 47979, 7, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 55096, 47214, 48844, 7, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 56041, 47979, 49779, 7, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 56986, 48844, 50824, 7, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 58141, 49779, 51946, 7, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 59296, 52408, 53584, 7, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 60472, 52996, 54340, 7, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 61648, 53584, 55096, 7, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 63160, 54340, 56041, 7, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 64672, 55096, 56986, 7, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 66562, 56041, 58141, 7, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 68452, 59296, 61648, 7, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 70412, 60472, 63160, 7, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 72372, 61648, 64672, 7, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 74892, 63160, 66562, 7, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 77412, 68452, 72372, 7, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 80352, 70412, 74892, 7, nmax);

        simdtrf::transform_i_inner(buffer, 83292, 80352, 15, 7, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 83292, 91, nmax);

        simdtrf::transform_i_inner(buffer, 83292, 77412, 15, 7, nmax);

        simdtrf::transform_g_outer(values + 819 * nvalues + n * npairs, nvalues, buffer, 83292,
                                   91, nmax);
    }

    for (size_t m = 0; m < 1638; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
