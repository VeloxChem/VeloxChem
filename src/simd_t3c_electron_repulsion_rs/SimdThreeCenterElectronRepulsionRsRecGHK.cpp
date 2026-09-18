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


#include "SimdThreeCenterElectronRepulsionRsRecGHK.hpp"

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
#include "SimdTransferDH.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferFH.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferGH.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_ghk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_ghk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 314467, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2970 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 314467, 241972, 18045, dimensions);

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
                                                            15, 16}, ncols, fj, i * nprim_b + j,
                                                            fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 23, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16}, ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 76, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 79, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 82, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 85, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 88, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 91, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 94, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 97, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 100, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 103, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 106, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 109, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 112, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 115, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 118, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 121, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 124, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 127, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 130, 0, 3, 7, 8,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 136, 0, 3, 8, 9,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 142, 0, 3, 9, 10,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 148, 0, 3, 10, 11,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 154, 0, 3, 11, 12,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 160, 0, 3, 12, 13,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 166, 0, 3, 13, 14,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 172, 0, 3, 14, 15,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 178, 0, 3, 15, 16,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 184, 0, 3, 16, 17,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 190, 0, 3, 17, 18,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 196, 0, 3, 18, 19,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 202, 0, 3, 19, 20,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 208, 0, 3, 20, 21,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 214, 0, 3, 24, 25,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 220, 0, 3, 25, 26,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 226, 0, 3, 26, 27,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 232, 0, 3, 27, 28,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 238, 0, 3, 28, 29,
                                                                       97, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 244, 0, 3, 29, 30,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 250, 0, 3, 30, 31,
                                                                       103, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 256, 0, 3, 31, 32,
                                                                       106, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 262, 0, 3, 32, 33,
                                                                       109, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 268, 0, 3, 33, 34,
                                                                       112, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 274, 0, 3, 34, 35,
                                                                       115, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 280, 0, 3, 35, 36,
                                                                       118, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 286, 0, 3, 36, 37,
                                                                       121, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 292, 0, 3, 37, 38,
                                                                       124, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 40, 43,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 308, 0, 3, 43, 46,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 46, 49,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 49, 52,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 52, 55,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 55, 58,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 58, 61,
                                                                       166, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 61, 64,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 64, 67,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 67, 70,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 70, 73,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 73, 76,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 76, 79,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 85, 88,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 88, 91,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 91, 94,
                                                                       226, 232, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 94, 97,
                                                                       232, 238, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 97,
                                                                       100, 238, 244, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 100,
                                                                       103, 244, 250, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 103,
                                                                       106, 250, 256, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 106,
                                                                       109, 256, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 109,
                                                                       112, 262, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 112,
                                                                       115, 268, 274, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 528, 0, 3, 115,
                                                                       118, 274, 280, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 538, 0, 3, 118,
                                                                       121, 280, 286, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 121,
                                                                       124, 286, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 558, 0, 3, 130,
                                                                       136, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 573, 0, 3, 136,
                                                                       142, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 142,
                                                                       148, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 603, 0, 3, 148,
                                                                       154, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 618, 0, 3, 154,
                                                                       160, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 633, 0, 3, 160,
                                                                       166, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 648, 0, 3, 166,
                                                                       172, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 663, 0, 3, 172,
                                                                       178, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 678, 0, 3, 178,
                                                                       184, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 693, 0, 3, 184,
                                                                       190, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 708, 0, 3, 190,
                                                                       196, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 723, 0, 3, 196,
                                                                       202, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 738, 0, 3, 214,
                                                                       220, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 753, 0, 3, 220,
                                                                       226, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 768, 0, 3, 226,
                                                                       232, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 783, 0, 3, 232,
                                                                       238, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 798, 0, 3, 238,
                                                                       244, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 813, 0, 3, 244,
                                                                       250, 478, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 828, 0, 3, 250,
                                                                       256, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 843, 0, 3, 256,
                                                                       262, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 858, 0, 3, 262,
                                                                       268, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 873, 0, 3, 268,
                                                                       274, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 888, 0, 3, 274,
                                                                       280, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 903, 0, 3, 280,
                                                                       286, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 918, 0, 3, 298,
                                                                       308, 558, 573, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 939, 0, 3, 308,
                                                                       318, 573, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 960, 0, 3, 318,
                                                                       328, 588, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 981, 0, 3, 328,
                                                                       338, 603, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 338,
                                                                       348, 618, 633, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1023, 0, 3, 348,
                                                                       358, 633, 648, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 358,
                                                                       368, 648, 663, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1065, 0, 3, 368,
                                                                       378, 663, 678, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 378,
                                                                       388, 678, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1107, 0, 3, 388,
                                                                       398, 693, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1128, 0, 3, 398,
                                                                       408, 708, 723, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1149, 0, 3, 428,
                                                                       438, 738, 753, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1170, 0, 3, 438,
                                                                       448, 753, 768, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1191, 0, 3, 448,
                                                                       458, 768, 783, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 458,
                                                                       468, 783, 798, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1233, 0, 3, 468,
                                                                       478, 798, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1254, 0, 3, 478,
                                                                       488, 813, 828, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1275, 0, 3, 488,
                                                                       498, 828, 843, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 498,
                                                                       508, 843, 858, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1317, 0, 3, 508,
                                                                       518, 858, 873, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1338, 0, 3, 518,
                                                                       528, 873, 888, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1359, 0, 3, 528,
                                                                       538, 888, 903, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 558,
                                                                       573, 918, 939, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 573,
                                                                       588, 939, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 588,
                                                                       603, 960, 981, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 603,
                                                                       618, 981, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 618,
                                                                       633, 1002, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 633,
                                                                       648, 1023, 1044, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 648,
                                                                       663, 1044, 1065, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 663,
                                                                       678, 1065, 1086, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 678,
                                                                       693, 1086, 1107, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 693,
                                                                       708, 1107, 1128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 738,
                                                                       753, 1149, 1170, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 753,
                                                                       768, 1170, 1191, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 768,
                                                                       783, 1191, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 783,
                                                                       798, 1212, 1233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 798,
                                                                       813, 1233, 1254, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 813,
                                                                       828, 1254, 1275, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 828,
                                                                       843, 1275, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 843,
                                                                       858, 1296, 1317, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 858,
                                                                       873, 1317, 1338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 873,
                                                                       888, 1338, 1359, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 918,
                                                                       939, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1976, 0, 3, 939,
                                                                       960, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2012, 0, 3, 960,
                                                                       981, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2048, 0, 3, 981,
                                                                       1002, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2084, 0, 3, 1002,
                                                                       1023, 1492, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2120, 0, 3, 1023,
                                                                       1044, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2156, 0, 3, 1044,
                                                                       1065, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 1065,
                                                                       1086, 1576, 1604, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2228, 0, 3, 1086,
                                                                       1107, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2264, 0, 3, 1149,
                                                                       1170, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2300, 0, 3, 1170,
                                                                       1191, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2336, 0, 3, 1191,
                                                                       1212, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1212,
                                                                       1233, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2408, 0, 3, 1233,
                                                                       1254, 1772, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2444, 0, 3, 1254,
                                                                       1275, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2480, 0, 3, 1275,
                                                                       1296, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2516, 0, 3, 1296,
                                                                       1317, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2552, 0, 3, 1317,
                                                                       1338, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2588, 0, 3, 1380,
                                                                       1408, 1940, 1976, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2633, 0, 3, 1408,
                                                                       1436, 1976, 2012, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2678, 0, 3, 1436,
                                                                       1464, 2012, 2048, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2723, 0, 3, 1464,
                                                                       1492, 2048, 2084, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2768, 0, 3, 1492,
                                                                       1520, 2084, 2120, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2813, 0, 3, 1520,
                                                                       1548, 2120, 2156, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2858, 0, 3, 1548,
                                                                       1576, 2156, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2903, 0, 3, 1576,
                                                                       1604, 2192, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2948, 0, 3, 1660,
                                                                       1688, 2264, 2300, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2993, 0, 3, 1688,
                                                                       1716, 2300, 2336, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3038, 0, 3, 1716,
                                                                       1744, 2336, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3083, 0, 3, 1744,
                                                                       1772, 2372, 2408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3128, 0, 3, 1772,
                                                                       1800, 2408, 2444, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3173, 0, 3, 1800,
                                                                       1828, 2444, 2480, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3218, 0, 3, 1828,
                                                                       1856, 2480, 2516, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3263, 0, 3, 1856,
                                                                       1884, 2516, 2552, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3308, 0, 3, 1940,
                                                                       1976, 2588, 2633, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3363, 0, 3, 1976,
                                                                       2012, 2633, 2678, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3418, 0, 3, 2012,
                                                                       2048, 2678, 2723, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3473, 0, 3, 2048,
                                                                       2084, 2723, 2768, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3528, 0, 3, 2084,
                                                                       2120, 2768, 2813, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3583, 0, 3, 2120,
                                                                       2156, 2813, 2858, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3638, 0, 3, 2156,
                                                                       2192, 2858, 2903, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3693, 0, 3, 2264,
                                                                       2300, 2948, 2993, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3748, 0, 3, 2300,
                                                                       2336, 2993, 3038, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3803, 0, 3, 2336,
                                                                       2372, 3038, 3083, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3858, 0, 3, 2372,
                                                                       2408, 3083, 3128, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3913, 0, 3, 2408,
                                                                       2444, 3128, 3173, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2444,
                                                                       2480, 3173, 3218, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4023, 0, 3, 2480,
                                                                       2516, 3218, 3263, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4078, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4081, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4084, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4087, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4090, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4093, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4096, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4099, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4102, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4105, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4108, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4111, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4114, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4117, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4120, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4123, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4126, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4129, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4132, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4135, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4138, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4141, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4144, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4147, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4150, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4153, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4156, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4159, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4162, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4165, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4168, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4171, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4174, 3, 9, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4183, 3, 10, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4192, 3, 11, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4201, 3, 12, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4210, 3, 13, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4219, 3, 14, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4228, 3, 15, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4237, 3, 16, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4246, 3, 17, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4255, 3, 18, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4264, 3, 19, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4273, 3, 20, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4282, 3, 21, 82,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4291, 3, 26, 91,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4300, 3, 27, 94,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4309, 3, 28, 97,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4318, 3, 29, 100,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4327, 3, 30, 103,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4336, 3, 31, 106,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4345, 3, 32, 109,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4354, 3, 33, 112,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4363, 3, 34, 115,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4372, 3, 35, 118,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4381, 3, 36, 121,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4390, 3, 37, 124,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4399, 3, 38, 127,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4408, 3, 40, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4426, 3, 43, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4444, 3, 46, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4462, 3, 49, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4480, 3, 52, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4498, 3, 55, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4516, 3, 58, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4534, 3, 61, 172,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4552, 3, 64, 178,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4570, 3, 67, 184,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4588, 3, 70, 190,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4606, 3, 73, 196,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4624, 3, 76, 202,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4642, 3, 79, 208,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4660, 3, 85, 214,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4678, 3, 88, 220,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4696, 3, 91, 226,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4714, 3, 94, 232,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4732, 3, 97, 238,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4750, 3, 100, 244,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4768, 3, 103, 250,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4786, 3, 106, 256,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4804, 3, 109, 262,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4822, 3, 112, 268,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4840, 3, 115, 274,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4858, 3, 118, 280,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4876, 3, 121, 286,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4894, 3, 124, 292,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4912, 3, 130, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4942, 3, 136, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4972, 3, 142, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5002, 3, 148, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5032, 3, 154, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5062, 3, 160, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5092, 3, 166, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5122, 3, 172, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5152, 3, 178, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5182, 3, 184, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5212, 3, 190, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5242, 3, 196, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5272, 3, 202, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5302, 3, 214, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5332, 3, 220, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5362, 3, 226, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5392, 3, 232, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5422, 3, 238, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5452, 3, 244, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5482, 3, 250, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5512, 3, 256, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5542, 3, 262, 508,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5572, 3, 268, 518,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5602, 3, 274, 528,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5632, 3, 280, 538,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5662, 3, 286, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5692, 3, 298, 558,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5737, 3, 308, 573,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5782, 3, 318, 588,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5827, 3, 328, 603,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5872, 3, 338, 618,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5917, 3, 348, 633,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5962, 3, 358, 648,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6007, 3, 368, 663,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6052, 3, 378, 678,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6097, 3, 388, 693,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6142, 3, 398, 708,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6187, 3, 408, 723,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6232, 3, 428, 738,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6277, 3, 438, 753,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6322, 3, 448, 768,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6367, 3, 458, 783,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6412, 3, 468, 798,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6457, 3, 478, 813,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6502, 3, 488, 828,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6547, 3, 498, 843,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6592, 3, 508, 858,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6637, 3, 518, 873,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6682, 3, 528, 888,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 6727, 3, 538, 903,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6772, 3, 558, 918,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6835, 3, 573, 939,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6898, 3, 588, 960,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6961, 3, 603, 981,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7024, 3, 618,
                                                                       1002, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7087, 3, 633,
                                                                       1023, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7150, 3, 648,
                                                                       1044, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7213, 3, 663,
                                                                       1065, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7276, 3, 678,
                                                                       1086, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7339, 3, 693,
                                                                       1107, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7402, 3, 708,
                                                                       1128, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7465, 3, 738,
                                                                       1149, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7528, 3, 753,
                                                                       1170, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7591, 3, 768,
                                                                       1191, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7654, 3, 783,
                                                                       1212, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7717, 3, 798,
                                                                       1233, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7780, 3, 813,
                                                                       1254, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7843, 3, 828,
                                                                       1275, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7906, 3, 843,
                                                                       1296, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 7969, 3, 858,
                                                                       1317, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8032, 3, 873,
                                                                       1338, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8095, 3, 888,
                                                                       1359, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8158, 3, 918,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8242, 3, 939,
                                                                       1408, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8326, 3, 960,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8410, 3, 981,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8494, 3, 1002,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8578, 3, 1023,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8662, 3, 1044,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8746, 3, 1065,
                                                                       1576, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8830, 3, 1086,
                                                                       1604, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8914, 3, 1107,
                                                                       1632, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 8998, 3, 1149,
                                                                       1660, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9082, 3, 1170,
                                                                       1688, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9166, 3, 1191,
                                                                       1716, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9250, 3, 1212,
                                                                       1744, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9334, 3, 1233,
                                                                       1772, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9418, 3, 1254,
                                                                       1800, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9502, 3, 1275,
                                                                       1828, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9586, 3, 1296,
                                                                       1856, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9670, 3, 1317,
                                                                       1884, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9754, 3, 1338,
                                                                       1912, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9838, 3, 1380,
                                                                       1940, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 9946, 3, 1408,
                                                                       1976, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10054, 3, 1436,
                                                                       2012, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10162, 3, 1464,
                                                                       2048, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10270, 3, 1492,
                                                                       2084, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10378, 3, 1520,
                                                                       2120, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10486, 3, 1548,
                                                                       2156, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10594, 3, 1576,
                                                                       2192, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10702, 3, 1604,
                                                                       2228, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10810, 3, 1660,
                                                                       2264, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 10918, 3, 1688,
                                                                       2300, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11026, 3, 1716,
                                                                       2336, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11134, 3, 1744,
                                                                       2372, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11242, 3, 1772,
                                                                       2408, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11350, 3, 1800,
                                                                       2444, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11458, 3, 1828,
                                                                       2480, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11566, 3, 1856,
                                                                       2516, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11674, 3, 1884,
                                                                       2552, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11782, 3, 1940,
                                                                       2588, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 11917, 3, 1976,
                                                                       2633, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12052, 3, 2012,
                                                                       2678, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12187, 3, 2048,
                                                                       2723, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12322, 3, 2084,
                                                                       2768, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12457, 3, 2120,
                                                                       2813, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12592, 3, 2156,
                                                                       2858, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12727, 3, 2192,
                                                                       2903, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12862, 3, 2264,
                                                                       2948, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 12997, 3, 2300,
                                                                       2993, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13132, 3, 2336,
                                                                       3038, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13267, 3, 2372,
                                                                       3083, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13402, 3, 2408,
                                                                       3128, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13537, 3, 2444,
                                                                       3173, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13672, 3, 2480,
                                                                       3218, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13807, 3, 2516,
                                                                       3263, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 13942, 3, 2588,
                                                                       3308, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14107, 3, 2633,
                                                                       3363, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14272, 3, 2678,
                                                                       3418, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14437, 3, 2723,
                                                                       3473, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14602, 3, 2768,
                                                                       3528, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14767, 3, 2813,
                                                                       3583, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 14932, 3, 2858,
                                                                       3638, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 15097, 3, 2948,
                                                                       3693, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 15262, 3, 2993,
                                                                       3748, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 15427, 3, 3038,
                                                                       3803, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 15592, 3, 3083,
                                                                       3858, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 15757, 3, 3128,
                                                                       3913, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 15922, 3, 3173,
                                                                       3968, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16087, 3, 3218,
                                                                       4023, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16252, 3, 7, 8,
                                                                       4084, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16258, 3, 8, 9,
                                                                       4087, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16264, 3, 9, 10,
                                                                       4090, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16270, 3, 10, 11,
                                                                       4093, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16276, 3, 11, 12,
                                                                       4096, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16282, 3, 12, 13,
                                                                       4099, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16288, 3, 13, 14,
                                                                       4102, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16294, 3, 14, 15,
                                                                       4105, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16300, 3, 15, 16,
                                                                       4108, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16306, 3, 16, 17,
                                                                       4111, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16312, 3, 17, 18,
                                                                       4114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16318, 3, 18, 19,
                                                                       4117, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16324, 3, 19, 20,
                                                                       4120, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16330, 3, 20, 21,
                                                                       4123, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16336, 3, 24, 25,
                                                                       4132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16342, 3, 25, 26,
                                                                       4135, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16348, 3, 26, 27,
                                                                       4138, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16354, 3, 27, 28,
                                                                       4141, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16360, 3, 28, 29,
                                                                       4144, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16366, 3, 29, 30,
                                                                       4147, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16372, 3, 30, 31,
                                                                       4150, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16378, 3, 31, 32,
                                                                       4153, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16384, 3, 32, 33,
                                                                       4156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16390, 3, 33, 34,
                                                                       4159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16396, 3, 34, 35,
                                                                       4162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16402, 3, 35, 36,
                                                                       4165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16408, 3, 36, 37,
                                                                       4168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 16414, 3, 37, 38,
                                                                       4171, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16420, 0, 3,
                                                                       16252, 4084, 16258, 4174,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16438, 0, 3,
                                                                       16258, 4087, 16264, 4183,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16456, 0, 3,
                                                                       16264, 4090, 16270, 4192,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16474, 0, 3,
                                                                       16270, 4093, 16276, 4201,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16492, 0, 3,
                                                                       16276, 4096, 16282, 4210,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16510, 0, 3,
                                                                       16282, 4099, 16288, 4219,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16528, 0, 3,
                                                                       16288, 4102, 16294, 4228,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16546, 0, 3,
                                                                       16294, 4105, 16300, 4237,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16564, 0, 3,
                                                                       16300, 4108, 16306, 4246,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16582, 0, 3,
                                                                       16306, 4111, 16312, 4255,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16600, 0, 3,
                                                                       16312, 4114, 16318, 4264,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16618, 0, 3,
                                                                       16318, 4117, 16324, 4273,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16636, 0, 3,
                                                                       16324, 4120, 16330, 4282,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16654, 0, 3,
                                                                       16336, 4132, 16342, 4291,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16672, 0, 3,
                                                                       16342, 4135, 16348, 4300,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16690, 0, 3,
                                                                       16348, 4138, 16354, 4309,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16708, 0, 3,
                                                                       16354, 4141, 16360, 4318,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16726, 0, 3,
                                                                       16360, 4144, 16366, 4327,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16744, 0, 3,
                                                                       16366, 4147, 16372, 4336,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16762, 0, 3,
                                                                       16372, 4150, 16378, 4345,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16780, 0, 3,
                                                                       16378, 4153, 16384, 4354,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16798, 0, 3,
                                                                       16384, 4156, 16390, 4363,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16816, 0, 3,
                                                                       16390, 4159, 16396, 4372,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16834, 0, 3,
                                                                       16396, 4162, 16402, 4381,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16852, 0, 3,
                                                                       16402, 4165, 16408, 4390,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 16870, 0, 3,
                                                                       16408, 4168, 16414, 4399,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 16888, 0, 3,
                                                                       16420, 4174, 16438, 130,
                                                                       136, 4444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 16924, 0, 3,
                                                                       16438, 4183, 16456, 136,
                                                                       142, 4462, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 16960, 0, 3,
                                                                       16456, 4192, 16474, 142,
                                                                       148, 4480, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 16996, 0, 3,
                                                                       16474, 4201, 16492, 148,
                                                                       154, 4498, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17032, 0, 3,
                                                                       16492, 4210, 16510, 154,
                                                                       160, 4516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17068, 0, 3,
                                                                       16510, 4219, 16528, 160,
                                                                       166, 4534, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17104, 0, 3,
                                                                       16528, 4228, 16546, 166,
                                                                       172, 4552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17140, 0, 3,
                                                                       16546, 4237, 16564, 172,
                                                                       178, 4570, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17176, 0, 3,
                                                                       16564, 4246, 16582, 178,
                                                                       184, 4588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17212, 0, 3,
                                                                       16582, 4255, 16600, 184,
                                                                       190, 4606, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17248, 0, 3,
                                                                       16600, 4264, 16618, 190,
                                                                       196, 4624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17284, 0, 3,
                                                                       16618, 4273, 16636, 196,
                                                                       202, 4642, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17320, 0, 3,
                                                                       16654, 4291, 16672, 214,
                                                                       220, 4696, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17356, 0, 3,
                                                                       16672, 4300, 16690, 220,
                                                                       226, 4714, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17392, 0, 3,
                                                                       16690, 4309, 16708, 226,
                                                                       232, 4732, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17428, 0, 3,
                                                                       16708, 4318, 16726, 232,
                                                                       238, 4750, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17464, 0, 3,
                                                                       16726, 4327, 16744, 238,
                                                                       244, 4768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17500, 0, 3,
                                                                       16744, 4336, 16762, 244,
                                                                       250, 4786, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17536, 0, 3,
                                                                       16762, 4345, 16780, 250,
                                                                       256, 4804, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17572, 0, 3,
                                                                       16780, 4354, 16798, 256,
                                                                       262, 4822, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17608, 0, 3,
                                                                       16798, 4363, 16816, 262,
                                                                       268, 4840, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17644, 0, 3,
                                                                       16816, 4372, 16834, 268,
                                                                       274, 4858, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17680, 0, 3,
                                                                       16834, 4381, 16852, 274,
                                                                       280, 4876, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 17716, 0, 3,
                                                                       16852, 4390, 16870, 280,
                                                                       286, 4894, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 17752, 0, 3,
                                                                       16888, 4444, 16924, 298,
                                                                       308, 4972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 17812, 0, 3,
                                                                       16924, 4462, 16960, 308,
                                                                       318, 5002, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 17872, 0, 3,
                                                                       16960, 4480, 16996, 318,
                                                                       328, 5032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 17932, 0, 3,
                                                                       16996, 4498, 17032, 328,
                                                                       338, 5062, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 17992, 0, 3,
                                                                       17032, 4516, 17068, 338,
                                                                       348, 5092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18052, 0, 3,
                                                                       17068, 4534, 17104, 348,
                                                                       358, 5122, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18112, 0, 3,
                                                                       17104, 4552, 17140, 358,
                                                                       368, 5152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18172, 0, 3,
                                                                       17140, 4570, 17176, 368,
                                                                       378, 5182, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18232, 0, 3,
                                                                       17176, 4588, 17212, 378,
                                                                       388, 5212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18292, 0, 3,
                                                                       17212, 4606, 17248, 388,
                                                                       398, 5242, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18352, 0, 3,
                                                                       17248, 4624, 17284, 398,
                                                                       408, 5272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18412, 0, 3,
                                                                       17320, 4696, 17356, 428,
                                                                       438, 5362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18472, 0, 3,
                                                                       17356, 4714, 17392, 438,
                                                                       448, 5392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18532, 0, 3,
                                                                       17392, 4732, 17428, 448,
                                                                       458, 5422, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18592, 0, 3,
                                                                       17428, 4750, 17464, 458,
                                                                       468, 5452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18652, 0, 3,
                                                                       17464, 4768, 17500, 468,
                                                                       478, 5482, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18712, 0, 3,
                                                                       17500, 4786, 17536, 478,
                                                                       488, 5512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18772, 0, 3,
                                                                       17536, 4804, 17572, 488,
                                                                       498, 5542, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18832, 0, 3,
                                                                       17572, 4822, 17608, 498,
                                                                       508, 5572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18892, 0, 3,
                                                                       17608, 4840, 17644, 508,
                                                                       518, 5602, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 18952, 0, 3,
                                                                       17644, 4858, 17680, 518,
                                                                       528, 5632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 19012, 0, 3,
                                                                       17680, 4876, 17716, 528,
                                                                       538, 5662, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19072, 0, 3,
                                                                       17752, 4972, 17812, 558,
                                                                       573, 5782, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19162, 0, 3,
                                                                       17812, 5002, 17872, 573,
                                                                       588, 5827, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19252, 0, 3,
                                                                       17872, 5032, 17932, 588,
                                                                       603, 5872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19342, 0, 3,
                                                                       17932, 5062, 17992, 603,
                                                                       618, 5917, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19432, 0, 3,
                                                                       17992, 5092, 18052, 618,
                                                                       633, 5962, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19522, 0, 3,
                                                                       18052, 5122, 18112, 633,
                                                                       648, 6007, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19612, 0, 3,
                                                                       18112, 5152, 18172, 648,
                                                                       663, 6052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19702, 0, 3,
                                                                       18172, 5182, 18232, 663,
                                                                       678, 6097, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19792, 0, 3,
                                                                       18232, 5212, 18292, 678,
                                                                       693, 6142, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19882, 0, 3,
                                                                       18292, 5242, 18352, 693,
                                                                       708, 6187, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 19972, 0, 3,
                                                                       18412, 5362, 18472, 738,
                                                                       753, 6322, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20062, 0, 3,
                                                                       18472, 5392, 18532, 753,
                                                                       768, 6367, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20152, 0, 3,
                                                                       18532, 5422, 18592, 768,
                                                                       783, 6412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20242, 0, 3,
                                                                       18592, 5452, 18652, 783,
                                                                       798, 6457, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20332, 0, 3,
                                                                       18652, 5482, 18712, 798,
                                                                       813, 6502, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20422, 0, 3,
                                                                       18712, 5512, 18772, 813,
                                                                       828, 6547, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20512, 0, 3,
                                                                       18772, 5542, 18832, 828,
                                                                       843, 6592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20602, 0, 3,
                                                                       18832, 5572, 18892, 843,
                                                                       858, 6637, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20692, 0, 3,
                                                                       18892, 5602, 18952, 858,
                                                                       873, 6682, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 20782, 0, 3,
                                                                       18952, 5632, 19012, 873,
                                                                       888, 6727, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 20872, 0, 3,
                                                                       19072, 5782, 19162, 918,
                                                                       939, 6898, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 20998, 0, 3,
                                                                       19162, 5827, 19252, 939,
                                                                       960, 6961, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21124, 0, 3,
                                                                       19252, 5872, 19342, 960,
                                                                       981, 7024, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21250, 0, 3,
                                                                       19342, 5917, 19432, 981,
                                                                       1002, 7087, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21376, 0, 3,
                                                                       19432, 5962, 19522, 1002,
                                                                       1023, 7150, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21502, 0, 3,
                                                                       19522, 6007, 19612, 1023,
                                                                       1044, 7213, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21628, 0, 3,
                                                                       19612, 6052, 19702, 1044,
                                                                       1065, 7276, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21754, 0, 3,
                                                                       19702, 6097, 19792, 1065,
                                                                       1086, 7339, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 21880, 0, 3,
                                                                       19792, 6142, 19882, 1086,
                                                                       1107, 7402, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22006, 0, 3,
                                                                       19972, 6322, 20062, 1149,
                                                                       1170, 7591, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22132, 0, 3,
                                                                       20062, 6367, 20152, 1170,
                                                                       1191, 7654, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22258, 0, 3,
                                                                       20152, 6412, 20242, 1191,
                                                                       1212, 7717, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22384, 0, 3,
                                                                       20242, 6457, 20332, 1212,
                                                                       1233, 7780, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22510, 0, 3,
                                                                       20332, 6502, 20422, 1233,
                                                                       1254, 7843, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22636, 0, 3,
                                                                       20422, 6547, 20512, 1254,
                                                                       1275, 7906, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22762, 0, 3,
                                                                       20512, 6592, 20602, 1275,
                                                                       1296, 7969, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 22888, 0, 3,
                                                                       20602, 6637, 20692, 1296,
                                                                       1317, 8032, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 23014, 0, 3,
                                                                       20692, 6682, 20782, 1317,
                                                                       1338, 8095, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23140, 0, 3,
                                                                       20872, 6898, 20998, 1380,
                                                                       1408, 8326, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23308, 0, 3,
                                                                       20998, 6961, 21124, 1408,
                                                                       1436, 8410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23476, 0, 3,
                                                                       21124, 7024, 21250, 1436,
                                                                       1464, 8494, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23644, 0, 3,
                                                                       21250, 7087, 21376, 1464,
                                                                       1492, 8578, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23812, 0, 3,
                                                                       21376, 7150, 21502, 1492,
                                                                       1520, 8662, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 23980, 0, 3,
                                                                       21502, 7213, 21628, 1520,
                                                                       1548, 8746, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24148, 0, 3,
                                                                       21628, 7276, 21754, 1548,
                                                                       1576, 8830, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24316, 0, 3,
                                                                       21754, 7339, 21880, 1576,
                                                                       1604, 8914, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24484, 0, 3,
                                                                       22006, 7591, 22132, 1660,
                                                                       1688, 9166, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24652, 0, 3,
                                                                       22132, 7654, 22258, 1688,
                                                                       1716, 9250, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24820, 0, 3,
                                                                       22258, 7717, 22384, 1716,
                                                                       1744, 9334, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 24988, 0, 3,
                                                                       22384, 7780, 22510, 1744,
                                                                       1772, 9418, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 25156, 0, 3,
                                                                       22510, 7843, 22636, 1772,
                                                                       1800, 9502, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 25324, 0, 3,
                                                                       22636, 7906, 22762, 1800,
                                                                       1828, 9586, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 25492, 0, 3,
                                                                       22762, 7969, 22888, 1828,
                                                                       1856, 9670, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 25660, 0, 3,
                                                                       22888, 8032, 23014, 1856,
                                                                       1884, 9754, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 25828, 0, 3,
                                                                       23140, 8326, 23308, 1940,
                                                                       1976, 10054, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26044, 0, 3,
                                                                       23308, 8410, 23476, 1976,
                                                                       2012, 10162, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26260, 0, 3,
                                                                       23476, 8494, 23644, 2012,
                                                                       2048, 10270, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26476, 0, 3,
                                                                       23644, 8578, 23812, 2048,
                                                                       2084, 10378, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26692, 0, 3,
                                                                       23812, 8662, 23980, 2084,
                                                                       2120, 10486, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 26908, 0, 3,
                                                                       23980, 8746, 24148, 2120,
                                                                       2156, 10594, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 27124, 0, 3,
                                                                       24148, 8830, 24316, 2156,
                                                                       2192, 10702, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 27340, 0, 3,
                                                                       24484, 9166, 24652, 2264,
                                                                       2300, 11026, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 27556, 0, 3,
                                                                       24652, 9250, 24820, 2300,
                                                                       2336, 11134, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 27772, 0, 3,
                                                                       24820, 9334, 24988, 2336,
                                                                       2372, 11242, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 27988, 0, 3,
                                                                       24988, 9418, 25156, 2372,
                                                                       2408, 11350, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 28204, 0, 3,
                                                                       25156, 9502, 25324, 2408,
                                                                       2444, 11458, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 28420, 0, 3,
                                                                       25324, 9586, 25492, 2444,
                                                                       2480, 11566, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 28636, 0, 3,
                                                                       25492, 9670, 25660, 2480,
                                                                       2516, 11674, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 28852, 0, 3,
                                                                       25828, 10054, 26044, 2588,
                                                                       2633, 12052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 29122, 0, 3,
                                                                       26044, 10162, 26260, 2633,
                                                                       2678, 12187, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 29392, 0, 3,
                                                                       26260, 10270, 26476, 2678,
                                                                       2723, 12322, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 29662, 0, 3,
                                                                       26476, 10378, 26692, 2723,
                                                                       2768, 12457, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 29932, 0, 3,
                                                                       26692, 10486, 26908, 2768,
                                                                       2813, 12592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 30202, 0, 3,
                                                                       26908, 10594, 27124, 2813,
                                                                       2858, 12727, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 30472, 0, 3,
                                                                       27340, 11026, 27556, 2948,
                                                                       2993, 13132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 30742, 0, 3,
                                                                       27556, 11134, 27772, 2993,
                                                                       3038, 13267, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 31012, 0, 3,
                                                                       27772, 11242, 27988, 3038,
                                                                       3083, 13402, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 31282, 0, 3,
                                                                       27988, 11350, 28204, 3083,
                                                                       3128, 13537, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 31552, 0, 3,
                                                                       28204, 11458, 28420, 3128,
                                                                       3173, 13672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 31822, 0, 3,
                                                                       28420, 11566, 28636, 3173,
                                                                       3218, 13807, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 32092, 0, 3,
                                                                       28852, 12052, 29122, 3308,
                                                                       3363, 14272, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 32422, 0, 3,
                                                                       29122, 12187, 29392, 3363,
                                                                       3418, 14437, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 32752, 0, 3,
                                                                       29392, 12322, 29662, 3418,
                                                                       3473, 14602, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 33082, 0, 3,
                                                                       29662, 12457, 29932, 3473,
                                                                       3528, 14767, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 33412, 0, 3,
                                                                       29932, 12592, 30202, 3528,
                                                                       3583, 14932, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 33742, 0, 3,
                                                                       30472, 13132, 30742, 3693,
                                                                       3748, 15427, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 34072, 0, 3,
                                                                       30742, 13267, 31012, 3748,
                                                                       3803, 15592, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 34402, 0, 3,
                                                                       31012, 13402, 31282, 3803,
                                                                       3858, 15757, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 34732, 0, 3,
                                                                       31282, 13537, 31552, 3858,
                                                                       3913, 15922, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 35062, 0, 3,
                                                                       31552, 13672, 31822, 3913,
                                                                       3968, 16087, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35392, 3, 4078,
                                                                       4081, 16252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35402, 3, 4081,
                                                                       4084, 16258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35412, 3, 4084,
                                                                       4087, 16264, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35422, 3, 4087,
                                                                       4090, 16270, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35432, 3, 4090,
                                                                       4093, 16276, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35442, 3, 4093,
                                                                       4096, 16282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35452, 3, 4096,
                                                                       4099, 16288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35462, 3, 4099,
                                                                       4102, 16294, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35472, 3, 4102,
                                                                       4105, 16300, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35482, 3, 4105,
                                                                       4108, 16306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35492, 3, 4108,
                                                                       4111, 16312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35502, 3, 4111,
                                                                       4114, 16318, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35512, 3, 4114,
                                                                       4117, 16324, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35522, 3, 4117,
                                                                       4120, 16330, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35532, 3, 4126,
                                                                       4129, 16336, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35542, 3, 4129,
                                                                       4132, 16342, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35552, 3, 4132,
                                                                       4135, 16348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35562, 3, 4135,
                                                                       4138, 16354, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35572, 3, 4138,
                                                                       4141, 16360, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35582, 3, 4141,
                                                                       4144, 16366, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35592, 3, 4144,
                                                                       4147, 16372, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35602, 3, 4147,
                                                                       4150, 16378, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35612, 3, 4150,
                                                                       4153, 16384, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35622, 3, 4153,
                                                                       4156, 16390, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35632, 3, 4156,
                                                                       4159, 16396, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35642, 3, 4159,
                                                                       4162, 16402, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35652, 3, 4162,
                                                                       4165, 16408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 35662, 3, 4165,
                                                                       4168, 16414, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35672, 0, 3,
                                                                       35392, 16252, 35402,
                                                                       16420, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35702, 0, 3,
                                                                       35402, 16258, 35412,
                                                                       16438, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35732, 0, 3,
                                                                       35412, 16264, 35422,
                                                                       16456, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35762, 0, 3,
                                                                       35422, 16270, 35432,
                                                                       16474, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35792, 0, 3,
                                                                       35432, 16276, 35442,
                                                                       16492, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35822, 0, 3,
                                                                       35442, 16282, 35452,
                                                                       16510, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35852, 0, 3,
                                                                       35452, 16288, 35462,
                                                                       16528, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35882, 0, 3,
                                                                       35462, 16294, 35472,
                                                                       16546, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35912, 0, 3,
                                                                       35472, 16300, 35482,
                                                                       16564, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35942, 0, 3,
                                                                       35482, 16306, 35492,
                                                                       16582, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 35972, 0, 3,
                                                                       35492, 16312, 35502,
                                                                       16600, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36002, 0, 3,
                                                                       35502, 16318, 35512,
                                                                       16618, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36032, 0, 3,
                                                                       35512, 16324, 35522,
                                                                       16636, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36062, 0, 3,
                                                                       35532, 16336, 35542,
                                                                       16654, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36092, 0, 3,
                                                                       35542, 16342, 35552,
                                                                       16672, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36122, 0, 3,
                                                                       35552, 16348, 35562,
                                                                       16690, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36152, 0, 3,
                                                                       35562, 16354, 35572,
                                                                       16708, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36182, 0, 3,
                                                                       35572, 16360, 35582,
                                                                       16726, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36212, 0, 3,
                                                                       35582, 16366, 35592,
                                                                       16744, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36242, 0, 3,
                                                                       35592, 16372, 35602,
                                                                       16762, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36272, 0, 3,
                                                                       35602, 16378, 35612,
                                                                       16780, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36302, 0, 3,
                                                                       35612, 16384, 35622,
                                                                       16798, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36332, 0, 3,
                                                                       35622, 16390, 35632,
                                                                       16816, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36362, 0, 3,
                                                                       35632, 16396, 35642,
                                                                       16834, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36392, 0, 3,
                                                                       35642, 16402, 35652,
                                                                       16852, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36422, 0, 3,
                                                                       35652, 16408, 35662,
                                                                       16870, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36452, 0, 3,
                                                                       35672, 16420, 35702, 4408,
                                                                       4426, 16888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36512, 0, 3,
                                                                       35702, 16438, 35732, 4426,
                                                                       4444, 16924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36572, 0, 3,
                                                                       35732, 16456, 35762, 4444,
                                                                       4462, 16960, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36632, 0, 3,
                                                                       35762, 16474, 35792, 4462,
                                                                       4480, 16996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36692, 0, 3,
                                                                       35792, 16492, 35822, 4480,
                                                                       4498, 17032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36752, 0, 3,
                                                                       35822, 16510, 35852, 4498,
                                                                       4516, 17068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36812, 0, 3,
                                                                       35852, 16528, 35882, 4516,
                                                                       4534, 17104, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36872, 0, 3,
                                                                       35882, 16546, 35912, 4534,
                                                                       4552, 17140, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36932, 0, 3,
                                                                       35912, 16564, 35942, 4552,
                                                                       4570, 17176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36992, 0, 3,
                                                                       35942, 16582, 35972, 4570,
                                                                       4588, 17212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37052, 0, 3,
                                                                       35972, 16600, 36002, 4588,
                                                                       4606, 17248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37112, 0, 3,
                                                                       36002, 16618, 36032, 4606,
                                                                       4624, 17284, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37172, 0, 3,
                                                                       36062, 16654, 36092, 4660,
                                                                       4678, 17320, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37232, 0, 3,
                                                                       36092, 16672, 36122, 4678,
                                                                       4696, 17356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37292, 0, 3,
                                                                       36122, 16690, 36152, 4696,
                                                                       4714, 17392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37352, 0, 3,
                                                                       36152, 16708, 36182, 4714,
                                                                       4732, 17428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37412, 0, 3,
                                                                       36182, 16726, 36212, 4732,
                                                                       4750, 17464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37472, 0, 3,
                                                                       36212, 16744, 36242, 4750,
                                                                       4768, 17500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37532, 0, 3,
                                                                       36242, 16762, 36272, 4768,
                                                                       4786, 17536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37592, 0, 3,
                                                                       36272, 16780, 36302, 4786,
                                                                       4804, 17572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37652, 0, 3,
                                                                       36302, 16798, 36332, 4804,
                                                                       4822, 17608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37712, 0, 3,
                                                                       36332, 16816, 36362, 4822,
                                                                       4840, 17644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37772, 0, 3,
                                                                       36362, 16834, 36392, 4840,
                                                                       4858, 17680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37832, 0, 3,
                                                                       36392, 16852, 36422, 4858,
                                                                       4876, 17716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 37892, 0, 3,
                                                                       36452, 16888, 36512, 4912,
                                                                       4942, 17752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 37992, 0, 3,
                                                                       36512, 16924, 36572, 4942,
                                                                       4972, 17812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38092, 0, 3,
                                                                       36572, 16960, 36632, 4972,
                                                                       5002, 17872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38192, 0, 3,
                                                                       36632, 16996, 36692, 5002,
                                                                       5032, 17932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38292, 0, 3,
                                                                       36692, 17032, 36752, 5032,
                                                                       5062, 17992, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38392, 0, 3,
                                                                       36752, 17068, 36812, 5062,
                                                                       5092, 18052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38492, 0, 3,
                                                                       36812, 17104, 36872, 5092,
                                                                       5122, 18112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38592, 0, 3,
                                                                       36872, 17140, 36932, 5122,
                                                                       5152, 18172, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38692, 0, 3,
                                                                       36932, 17176, 36992, 5152,
                                                                       5182, 18232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38792, 0, 3,
                                                                       36992, 17212, 37052, 5182,
                                                                       5212, 18292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38892, 0, 3,
                                                                       37052, 17248, 37112, 5212,
                                                                       5242, 18352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38992, 0, 3,
                                                                       37172, 17320, 37232, 5302,
                                                                       5332, 18412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39092, 0, 3,
                                                                       37232, 17356, 37292, 5332,
                                                                       5362, 18472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39192, 0, 3,
                                                                       37292, 17392, 37352, 5362,
                                                                       5392, 18532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39292, 0, 3,
                                                                       37352, 17428, 37412, 5392,
                                                                       5422, 18592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39392, 0, 3,
                                                                       37412, 17464, 37472, 5422,
                                                                       5452, 18652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39492, 0, 3,
                                                                       37472, 17500, 37532, 5452,
                                                                       5482, 18712, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39592, 0, 3,
                                                                       37532, 17536, 37592, 5482,
                                                                       5512, 18772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39692, 0, 3,
                                                                       37592, 17572, 37652, 5512,
                                                                       5542, 18832, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39792, 0, 3,
                                                                       37652, 17608, 37712, 5542,
                                                                       5572, 18892, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39892, 0, 3,
                                                                       37712, 17644, 37772, 5572,
                                                                       5602, 18952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 39992, 0, 3,
                                                                       37772, 17680, 37832, 5602,
                                                                       5632, 19012, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40092, 0, 3,
                                                                       37892, 17752, 37992, 5692,
                                                                       5737, 19072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40242, 0, 3,
                                                                       37992, 17812, 38092, 5737,
                                                                       5782, 19162, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40392, 0, 3,
                                                                       38092, 17872, 38192, 5782,
                                                                       5827, 19252, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40542, 0, 3,
                                                                       38192, 17932, 38292, 5827,
                                                                       5872, 19342, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40692, 0, 3,
                                                                       38292, 17992, 38392, 5872,
                                                                       5917, 19432, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40842, 0, 3,
                                                                       38392, 18052, 38492, 5917,
                                                                       5962, 19522, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40992, 0, 3,
                                                                       38492, 18112, 38592, 5962,
                                                                       6007, 19612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41142, 0, 3,
                                                                       38592, 18172, 38692, 6007,
                                                                       6052, 19702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41292, 0, 3,
                                                                       38692, 18232, 38792, 6052,
                                                                       6097, 19792, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41442, 0, 3,
                                                                       38792, 18292, 38892, 6097,
                                                                       6142, 19882, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41592, 0, 3,
                                                                       38992, 18412, 39092, 6232,
                                                                       6277, 19972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41742, 0, 3,
                                                                       39092, 18472, 39192, 6277,
                                                                       6322, 20062, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 41892, 0, 3,
                                                                       39192, 18532, 39292, 6322,
                                                                       6367, 20152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42042, 0, 3,
                                                                       39292, 18592, 39392, 6367,
                                                                       6412, 20242, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42192, 0, 3,
                                                                       39392, 18652, 39492, 6412,
                                                                       6457, 20332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42342, 0, 3,
                                                                       39492, 18712, 39592, 6457,
                                                                       6502, 20422, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42492, 0, 3,
                                                                       39592, 18772, 39692, 6502,
                                                                       6547, 20512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42642, 0, 3,
                                                                       39692, 18832, 39792, 6547,
                                                                       6592, 20602, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42792, 0, 3,
                                                                       39792, 18892, 39892, 6592,
                                                                       6637, 20692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 42942, 0, 3,
                                                                       39892, 18952, 39992, 6637,
                                                                       6682, 20782, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 43092, 0, 3,
                                                                       40092, 19072, 40242, 6772,
                                                                       6835, 20872, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 43302, 0, 3,
                                                                       40242, 19162, 40392, 6835,
                                                                       6898, 20998, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 43512, 0, 3,
                                                                       40392, 19252, 40542, 6898,
                                                                       6961, 21124, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 43722, 0, 3,
                                                                       40542, 19342, 40692, 6961,
                                                                       7024, 21250, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 43932, 0, 3,
                                                                       40692, 19432, 40842, 7024,
                                                                       7087, 21376, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 44142, 0, 3,
                                                                       40842, 19522, 40992, 7087,
                                                                       7150, 21502, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 44352, 0, 3,
                                                                       40992, 19612, 41142, 7150,
                                                                       7213, 21628, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 44562, 0, 3,
                                                                       41142, 19702, 41292, 7213,
                                                                       7276, 21754, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 44772, 0, 3,
                                                                       41292, 19792, 41442, 7276,
                                                                       7339, 21880, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 44982, 0, 3,
                                                                       41592, 19972, 41742, 7465,
                                                                       7528, 22006, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 45192, 0, 3,
                                                                       41742, 20062, 41892, 7528,
                                                                       7591, 22132, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 45402, 0, 3,
                                                                       41892, 20152, 42042, 7591,
                                                                       7654, 22258, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 45612, 0, 3,
                                                                       42042, 20242, 42192, 7654,
                                                                       7717, 22384, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 45822, 0, 3,
                                                                       42192, 20332, 42342, 7717,
                                                                       7780, 22510, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 46032, 0, 3,
                                                                       42342, 20422, 42492, 7780,
                                                                       7843, 22636, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 46242, 0, 3,
                                                                       42492, 20512, 42642, 7843,
                                                                       7906, 22762, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 46452, 0, 3,
                                                                       42642, 20602, 42792, 7906,
                                                                       7969, 22888, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 46662, 0, 3,
                                                                       42792, 20692, 42942, 7969,
                                                                       8032, 23014, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 46872, 0, 3,
                                                                       43092, 20872, 43302, 8158,
                                                                       8242, 23140, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 47152, 0, 3,
                                                                       43302, 20998, 43512, 8242,
                                                                       8326, 23308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 47432, 0, 3,
                                                                       43512, 21124, 43722, 8326,
                                                                       8410, 23476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 47712, 0, 3,
                                                                       43722, 21250, 43932, 8410,
                                                                       8494, 23644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 47992, 0, 3,
                                                                       43932, 21376, 44142, 8494,
                                                                       8578, 23812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 48272, 0, 3,
                                                                       44142, 21502, 44352, 8578,
                                                                       8662, 23980, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 48552, 0, 3,
                                                                       44352, 21628, 44562, 8662,
                                                                       8746, 24148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 48832, 0, 3,
                                                                       44562, 21754, 44772, 8746,
                                                                       8830, 24316, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 49112, 0, 3,
                                                                       44982, 22006, 45192, 8998,
                                                                       9082, 24484, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 49392, 0, 3,
                                                                       45192, 22132, 45402, 9082,
                                                                       9166, 24652, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 49672, 0, 3,
                                                                       45402, 22258, 45612, 9166,
                                                                       9250, 24820, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 49952, 0, 3,
                                                                       45612, 22384, 45822, 9250,
                                                                       9334, 24988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 50232, 0, 3,
                                                                       45822, 22510, 46032, 9334,
                                                                       9418, 25156, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 50512, 0, 3,
                                                                       46032, 22636, 46242, 9418,
                                                                       9502, 25324, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 50792, 0, 3,
                                                                       46242, 22762, 46452, 9502,
                                                                       9586, 25492, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 51072, 0, 3,
                                                                       46452, 22888, 46662, 9586,
                                                                       9670, 25660, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 51352, 0, 3,
                                                                       46872, 23140, 47152, 9838,
                                                                       9946, 25828, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 51712, 0, 3,
                                                                       47152, 23308, 47432, 9946,
                                                                       10054, 26044, ncols,
                                                                       gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 52072, 0, 3,
                                                                       47432, 23476, 47712,
                                                                       10054, 10162, 26260,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 52432, 0, 3,
                                                                       47712, 23644, 47992,
                                                                       10162, 10270, 26476,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 52792, 0, 3,
                                                                       47992, 23812, 48272,
                                                                       10270, 10378, 26692,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 53152, 0, 3,
                                                                       48272, 23980, 48552,
                                                                       10378, 10486, 26908,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 53512, 0, 3,
                                                                       48552, 24148, 48832,
                                                                       10486, 10594, 27124,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 53872, 0, 3,
                                                                       49112, 24484, 49392,
                                                                       10810, 10918, 27340,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 54232, 0, 3,
                                                                       49392, 24652, 49672,
                                                                       10918, 11026, 27556,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 54592, 0, 3,
                                                                       49672, 24820, 49952,
                                                                       11026, 11134, 27772,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 54952, 0, 3,
                                                                       49952, 24988, 50232,
                                                                       11134, 11242, 27988,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 55312, 0, 3,
                                                                       50232, 25156, 50512,
                                                                       11242, 11350, 28204,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 55672, 0, 3,
                                                                       50512, 25324, 50792,
                                                                       11350, 11458, 28420,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 56032, 0, 3,
                                                                       50792, 25492, 51072,
                                                                       11458, 11566, 28636,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 56392, 0, 3,
                                                                       51352, 25828, 51712,
                                                                       11782, 11917, 28852,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 56842, 0, 3,
                                                                       51712, 26044, 52072,
                                                                       11917, 12052, 29122,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 57292, 0, 3,
                                                                       52072, 26260, 52432,
                                                                       12052, 12187, 29392,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 57742, 0, 3,
                                                                       52432, 26476, 52792,
                                                                       12187, 12322, 29662,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 58192, 0, 3,
                                                                       52792, 26692, 53152,
                                                                       12322, 12457, 29932,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 58642, 0, 3,
                                                                       53152, 26908, 53512,
                                                                       12457, 12592, 30202,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 59092, 0, 3,
                                                                       53872, 27340, 54232,
                                                                       12862, 12997, 30472,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 59542, 0, 3,
                                                                       54232, 27556, 54592,
                                                                       12997, 13132, 30742,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 59992, 0, 3,
                                                                       54592, 27772, 54952,
                                                                       13132, 13267, 31012,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 60442, 0, 3,
                                                                       54952, 27988, 55312,
                                                                       13267, 13402, 31282,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 60892, 0, 3,
                                                                       55312, 28204, 55672,
                                                                       13402, 13537, 31552,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 61342, 0, 3,
                                                                       55672, 28420, 56032,
                                                                       13537, 13672, 31822,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 61792, 0, 3,
                                                                       56392, 28852, 56842,
                                                                       13942, 14107, 32092,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 62342, 0, 3,
                                                                       56842, 29122, 57292,
                                                                       14107, 14272, 32422,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 62892, 0, 3,
                                                                       57292, 29392, 57742,
                                                                       14272, 14437, 32752,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 63442, 0, 3,
                                                                       57742, 29662, 58192,
                                                                       14437, 14602, 33082,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 63992, 0, 3,
                                                                       58192, 29932, 58642,
                                                                       14602, 14767, 33412,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 64542, 0, 3,
                                                                       59092, 30472, 59542,
                                                                       15097, 15262, 33742,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 65092, 0, 3,
                                                                       59542, 30742, 59992,
                                                                       15262, 15427, 34072,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 65642, 0, 3,
                                                                       59992, 31012, 60442,
                                                                       15427, 15592, 34402,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 66192, 0, 3,
                                                                       60442, 31282, 60892,
                                                                       15592, 15757, 34732,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 66742, 0, 3,
                                                                       60892, 31552, 61342,
                                                                       15757, 15922, 35062,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67292, 3, 16252,
                                                                       16258, 35412, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67307, 3, 16258,
                                                                       16264, 35422, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67322, 3, 16264,
                                                                       16270, 35432, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67337, 3, 16270,
                                                                       16276, 35442, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67352, 3, 16276,
                                                                       16282, 35452, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67367, 3, 16282,
                                                                       16288, 35462, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67382, 3, 16288,
                                                                       16294, 35472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67397, 3, 16294,
                                                                       16300, 35482, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67412, 3, 16300,
                                                                       16306, 35492, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67427, 3, 16306,
                                                                       16312, 35502, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67442, 3, 16312,
                                                                       16318, 35512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67457, 3, 16318,
                                                                       16324, 35522, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67472, 3, 16336,
                                                                       16342, 35552, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67487, 3, 16342,
                                                                       16348, 35562, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67502, 3, 16348,
                                                                       16354, 35572, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67517, 3, 16354,
                                                                       16360, 35582, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67532, 3, 16360,
                                                                       16366, 35592, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67547, 3, 16366,
                                                                       16372, 35602, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67562, 3, 16372,
                                                                       16378, 35612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67577, 3, 16378,
                                                                       16384, 35622, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67592, 3, 16384,
                                                                       16390, 35632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67607, 3, 16390,
                                                                       16396, 35642, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67622, 3, 16396,
                                                                       16402, 35652, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 67637, 3, 16402,
                                                                       16408, 35662, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67652, 0, 3,
                                                                       67292, 35412, 67307,
                                                                       16420, 16438, 35732,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67697, 0, 3,
                                                                       67307, 35422, 67322,
                                                                       16438, 16456, 35762,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67742, 0, 3,
                                                                       67322, 35432, 67337,
                                                                       16456, 16474, 35792,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67787, 0, 3,
                                                                       67337, 35442, 67352,
                                                                       16474, 16492, 35822,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67832, 0, 3,
                                                                       67352, 35452, 67367,
                                                                       16492, 16510, 35852,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67877, 0, 3,
                                                                       67367, 35462, 67382,
                                                                       16510, 16528, 35882,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67922, 0, 3,
                                                                       67382, 35472, 67397,
                                                                       16528, 16546, 35912,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 67967, 0, 3,
                                                                       67397, 35482, 67412,
                                                                       16546, 16564, 35942,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68012, 0, 3,
                                                                       67412, 35492, 67427,
                                                                       16564, 16582, 35972,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68057, 0, 3,
                                                                       67427, 35502, 67442,
                                                                       16582, 16600, 36002,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68102, 0, 3,
                                                                       67442, 35512, 67457,
                                                                       16600, 16618, 36032,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68147, 0, 3,
                                                                       67472, 35552, 67487,
                                                                       16654, 16672, 36122,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68192, 0, 3,
                                                                       67487, 35562, 67502,
                                                                       16672, 16690, 36152,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68237, 0, 3,
                                                                       67502, 35572, 67517,
                                                                       16690, 16708, 36182,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68282, 0, 3,
                                                                       67517, 35582, 67532,
                                                                       16708, 16726, 36212,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68327, 0, 3,
                                                                       67532, 35592, 67547,
                                                                       16726, 16744, 36242,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68372, 0, 3,
                                                                       67547, 35602, 67562,
                                                                       16744, 16762, 36272,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68417, 0, 3,
                                                                       67562, 35612, 67577,
                                                                       16762, 16780, 36302,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68462, 0, 3,
                                                                       67577, 35622, 67592,
                                                                       16780, 16798, 36332,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68507, 0, 3,
                                                                       67592, 35632, 67607,
                                                                       16798, 16816, 36362,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68552, 0, 3,
                                                                       67607, 35642, 67622,
                                                                       16816, 16834, 36392,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 68597, 0, 3,
                                                                       67622, 35652, 67637,
                                                                       16834, 16852, 36422,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 68642, 0, 3,
                                                                       67652, 35732, 67697,
                                                                       16888, 16924, 36572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 68732, 0, 3,
                                                                       67697, 35762, 67742,
                                                                       16924, 16960, 36632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 68822, 0, 3,
                                                                       67742, 35792, 67787,
                                                                       16960, 16996, 36692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 68912, 0, 3,
                                                                       67787, 35822, 67832,
                                                                       16996, 17032, 36752,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69002, 0, 3,
                                                                       67832, 35852, 67877,
                                                                       17032, 17068, 36812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69092, 0, 3,
                                                                       67877, 35882, 67922,
                                                                       17068, 17104, 36872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69182, 0, 3,
                                                                       67922, 35912, 67967,
                                                                       17104, 17140, 36932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69272, 0, 3,
                                                                       67967, 35942, 68012,
                                                                       17140, 17176, 36992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69362, 0, 3,
                                                                       68012, 35972, 68057,
                                                                       17176, 17212, 37052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69452, 0, 3,
                                                                       68057, 36002, 68102,
                                                                       17212, 17248, 37112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69542, 0, 3,
                                                                       68147, 36122, 68192,
                                                                       17320, 17356, 37292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69632, 0, 3,
                                                                       68192, 36152, 68237,
                                                                       17356, 17392, 37352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69722, 0, 3,
                                                                       68237, 36182, 68282,
                                                                       17392, 17428, 37412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69812, 0, 3,
                                                                       68282, 36212, 68327,
                                                                       17428, 17464, 37472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69902, 0, 3,
                                                                       68327, 36242, 68372,
                                                                       17464, 17500, 37532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 69992, 0, 3,
                                                                       68372, 36272, 68417,
                                                                       17500, 17536, 37592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 70082, 0, 3,
                                                                       68417, 36302, 68462,
                                                                       17536, 17572, 37652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 70172, 0, 3,
                                                                       68462, 36332, 68507,
                                                                       17572, 17608, 37712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 70262, 0, 3,
                                                                       68507, 36362, 68552,
                                                                       17608, 17644, 37772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 70352, 0, 3,
                                                                       68552, 36392, 68597,
                                                                       17644, 17680, 37832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 70442, 0, 3,
                                                                       68642, 36572, 68732,
                                                                       17752, 17812, 38092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 70592, 0, 3,
                                                                       68732, 36632, 68822,
                                                                       17812, 17872, 38192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 70742, 0, 3,
                                                                       68822, 36692, 68912,
                                                                       17872, 17932, 38292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 70892, 0, 3,
                                                                       68912, 36752, 69002,
                                                                       17932, 17992, 38392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71042, 0, 3,
                                                                       69002, 36812, 69092,
                                                                       17992, 18052, 38492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71192, 0, 3,
                                                                       69092, 36872, 69182,
                                                                       18052, 18112, 38592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71342, 0, 3,
                                                                       69182, 36932, 69272,
                                                                       18112, 18172, 38692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71492, 0, 3,
                                                                       69272, 36992, 69362,
                                                                       18172, 18232, 38792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71642, 0, 3,
                                                                       69362, 37052, 69452,
                                                                       18232, 18292, 38892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71792, 0, 3,
                                                                       69542, 37292, 69632,
                                                                       18412, 18472, 39192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 71942, 0, 3,
                                                                       69632, 37352, 69722,
                                                                       18472, 18532, 39292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72092, 0, 3,
                                                                       69722, 37412, 69812,
                                                                       18532, 18592, 39392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72242, 0, 3,
                                                                       69812, 37472, 69902,
                                                                       18592, 18652, 39492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72392, 0, 3,
                                                                       69902, 37532, 69992,
                                                                       18652, 18712, 39592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72542, 0, 3,
                                                                       69992, 37592, 70082,
                                                                       18712, 18772, 39692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72692, 0, 3,
                                                                       70082, 37652, 70172,
                                                                       18772, 18832, 39792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72842, 0, 3,
                                                                       70172, 37712, 70262,
                                                                       18832, 18892, 39892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 72992, 0, 3,
                                                                       70262, 37772, 70352,
                                                                       18892, 18952, 39992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 73142, 0, 3,
                                                                       70442, 38092, 70592,
                                                                       19072, 19162, 40392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 73367, 0, 3,
                                                                       70592, 38192, 70742,
                                                                       19162, 19252, 40542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 73592, 0, 3,
                                                                       70742, 38292, 70892,
                                                                       19252, 19342, 40692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 73817, 0, 3,
                                                                       70892, 38392, 71042,
                                                                       19342, 19432, 40842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 74042, 0, 3,
                                                                       71042, 38492, 71192,
                                                                       19432, 19522, 40992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 74267, 0, 3,
                                                                       71192, 38592, 71342,
                                                                       19522, 19612, 41142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 74492, 0, 3,
                                                                       71342, 38692, 71492,
                                                                       19612, 19702, 41292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 74717, 0, 3,
                                                                       71492, 38792, 71642,
                                                                       19702, 19792, 41442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 74942, 0, 3,
                                                                       71792, 39192, 71942,
                                                                       19972, 20062, 41892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 75167, 0, 3,
                                                                       71942, 39292, 72092,
                                                                       20062, 20152, 42042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 75392, 0, 3,
                                                                       72092, 39392, 72242,
                                                                       20152, 20242, 42192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 75617, 0, 3,
                                                                       72242, 39492, 72392,
                                                                       20242, 20332, 42342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 75842, 0, 3,
                                                                       72392, 39592, 72542,
                                                                       20332, 20422, 42492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 76067, 0, 3,
                                                                       72542, 39692, 72692,
                                                                       20422, 20512, 42642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 76292, 0, 3,
                                                                       72692, 39792, 72842,
                                                                       20512, 20602, 42792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 76517, 0, 3,
                                                                       72842, 39892, 72992,
                                                                       20602, 20692, 42942,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 76742, 0, 3,
                                                                       73142, 40392, 73367,
                                                                       20872, 20998, 43512,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 77057, 0, 3,
                                                                       73367, 40542, 73592,
                                                                       20998, 21124, 43722,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 77372, 0, 3,
                                                                       73592, 40692, 73817,
                                                                       21124, 21250, 43932,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 77687, 0, 3,
                                                                       73817, 40842, 74042,
                                                                       21250, 21376, 44142,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 78002, 0, 3,
                                                                       74042, 40992, 74267,
                                                                       21376, 21502, 44352,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 78317, 0, 3,
                                                                       74267, 41142, 74492,
                                                                       21502, 21628, 44562,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 78632, 0, 3,
                                                                       74492, 41292, 74717,
                                                                       21628, 21754, 44772,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 78947, 0, 3,
                                                                       74942, 41892, 75167,
                                                                       22006, 22132, 45402,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 79262, 0, 3,
                                                                       75167, 42042, 75392,
                                                                       22132, 22258, 45612,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 79577, 0, 3,
                                                                       75392, 42192, 75617,
                                                                       22258, 22384, 45822,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 79892, 0, 3,
                                                                       75617, 42342, 75842,
                                                                       22384, 22510, 46032,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 80207, 0, 3,
                                                                       75842, 42492, 76067,
                                                                       22510, 22636, 46242,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 80522, 0, 3,
                                                                       76067, 42642, 76292,
                                                                       22636, 22762, 46452,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 80837, 0, 3,
                                                                       76292, 42792, 76517,
                                                                       22762, 22888, 46662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 81152, 0, 3,
                                                                       76742, 43512, 77057,
                                                                       23140, 23308, 47432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 81572, 0, 3,
                                                                       77057, 43722, 77372,
                                                                       23308, 23476, 47712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 81992, 0, 3,
                                                                       77372, 43932, 77687,
                                                                       23476, 23644, 47992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 82412, 0, 3,
                                                                       77687, 44142, 78002,
                                                                       23644, 23812, 48272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 82832, 0, 3,
                                                                       78002, 44352, 78317,
                                                                       23812, 23980, 48552,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 83252, 0, 3,
                                                                       78317, 44562, 78632,
                                                                       23980, 24148, 48832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 83672, 0, 3,
                                                                       78947, 45402, 79262,
                                                                       24484, 24652, 49672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 84092, 0, 3,
                                                                       79262, 45612, 79577,
                                                                       24652, 24820, 49952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 84512, 0, 3,
                                                                       79577, 45822, 79892,
                                                                       24820, 24988, 50232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 84932, 0, 3,
                                                                       79892, 46032, 80207,
                                                                       24988, 25156, 50512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 85352, 0, 3,
                                                                       80207, 46242, 80522,
                                                                       25156, 25324, 50792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 85772, 0, 3,
                                                                       80522, 46452, 80837,
                                                                       25324, 25492, 51072,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 86192, 0, 3,
                                                                       81152, 47432, 81572,
                                                                       25828, 26044, 52072,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 86732, 0, 3,
                                                                       81572, 47712, 81992,
                                                                       26044, 26260, 52432,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 87272, 0, 3,
                                                                       81992, 47992, 82412,
                                                                       26260, 26476, 52792,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 87812, 0, 3,
                                                                       82412, 48272, 82832,
                                                                       26476, 26692, 53152,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 88352, 0, 3,
                                                                       82832, 48552, 83252,
                                                                       26692, 26908, 53512,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 88892, 0, 3,
                                                                       83672, 49672, 84092,
                                                                       27340, 27556, 54592,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 89432, 0, 3,
                                                                       84092, 49952, 84512,
                                                                       27556, 27772, 54952,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 89972, 0, 3,
                                                                       84512, 50232, 84932,
                                                                       27772, 27988, 55312,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 90512, 0, 3,
                                                                       84932, 50512, 85352,
                                                                       27988, 28204, 55672,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 91052, 0, 3,
                                                                       85352, 50792, 85772,
                                                                       28204, 28420, 56032,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 91592, 0, 3,
                                                                       86192, 52072, 86732,
                                                                       28852, 29122, 57292,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 92267, 0, 3,
                                                                       86732, 52432, 87272,
                                                                       29122, 29392, 57742,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 92942, 0, 3,
                                                                       87272, 52792, 87812,
                                                                       29392, 29662, 58192,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 93617, 0, 3,
                                                                       87812, 53152, 88352,
                                                                       29662, 29932, 58642,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 94292, 0, 3,
                                                                       88892, 54592, 89432,
                                                                       30472, 30742, 59992,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 94967, 0, 3,
                                                                       89432, 54952, 89972,
                                                                       30742, 31012, 60442,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 95642, 0, 3,
                                                                       89972, 55312, 90512,
                                                                       31012, 31282, 60892,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 96317, 0, 3,
                                                                       90512, 55672, 91052,
                                                                       31282, 31552, 61342,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 96992, 0, 3,
                                                                       91592, 57292, 92267,
                                                                       32092, 32422, 62892,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 97817, 0, 3,
                                                                       92267, 57742, 92942,
                                                                       32422, 32752, 63442,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 98642, 0, 3,
                                                                       92942, 58192, 93617,
                                                                       32752, 33082, 63992,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 99467, 0, 3,
                                                                       94292, 59992, 94967,
                                                                       33742, 34072, 65642,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 100292, 0, 3,
                                                                       94967, 60442, 95642,
                                                                       34072, 34402, 66192,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 101117, 0, 3,
                                                                       95642, 60892, 96317,
                                                                       34402, 34732, 66742,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 101942, 3, 35392,
                                                                       35402, 67292, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 101963, 3, 35402,
                                                                       35412, 67307, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 101984, 3, 35412,
                                                                       35422, 67322, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102005, 3, 35422,
                                                                       35432, 67337, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102026, 3, 35432,
                                                                       35442, 67352, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102047, 3, 35442,
                                                                       35452, 67367, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102068, 3, 35452,
                                                                       35462, 67382, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102089, 3, 35462,
                                                                       35472, 67397, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102110, 3, 35472,
                                                                       35482, 67412, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102131, 3, 35482,
                                                                       35492, 67427, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102152, 3, 35492,
                                                                       35502, 67442, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102173, 3, 35502,
                                                                       35512, 67457, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102194, 3, 35532,
                                                                       35542, 67472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102215, 3, 35542,
                                                                       35552, 67487, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102236, 3, 35552,
                                                                       35562, 67502, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102257, 3, 35562,
                                                                       35572, 67517, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102278, 3, 35572,
                                                                       35582, 67532, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102299, 3, 35582,
                                                                       35592, 67547, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102320, 3, 35592,
                                                                       35602, 67562, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102341, 3, 35602,
                                                                       35612, 67577, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102362, 3, 35612,
                                                                       35622, 67592, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102383, 3, 35622,
                                                                       35632, 67607, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102404, 3, 35632,
                                                                       35642, 67622, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 102425, 3, 35642,
                                                                       35652, 67637, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 102446, 0, 3,
                                                                       101942, 67292, 101963,
                                                                       35672, 35702, 67652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 102509, 0, 3,
                                                                       101963, 67307, 101984,
                                                                       35702, 35732, 67697,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 102572, 0, 3,
                                                                       101984, 67322, 102005,
                                                                       35732, 35762, 67742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 102635, 0, 3,
                                                                       102005, 67337, 102026,
                                                                       35762, 35792, 67787,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 102698, 0, 3,
                                                                       102026, 67352, 102047,
                                                                       35792, 35822, 67832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 102761, 0, 3,
                                                                       102047, 67367, 102068,
                                                                       35822, 35852, 67877,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 102824, 0, 3,
                                                                       102068, 67382, 102089,
                                                                       35852, 35882, 67922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 102887, 0, 3,
                                                                       102089, 67397, 102110,
                                                                       35882, 35912, 67967,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 102950, 0, 3,
                                                                       102110, 67412, 102131,
                                                                       35912, 35942, 68012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103013, 0, 3,
                                                                       102131, 67427, 102152,
                                                                       35942, 35972, 68057,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103076, 0, 3,
                                                                       102152, 67442, 102173,
                                                                       35972, 36002, 68102,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103139, 0, 3,
                                                                       102194, 67472, 102215,
                                                                       36062, 36092, 68147,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103202, 0, 3,
                                                                       102215, 67487, 102236,
                                                                       36092, 36122, 68192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103265, 0, 3,
                                                                       102236, 67502, 102257,
                                                                       36122, 36152, 68237,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103328, 0, 3,
                                                                       102257, 67517, 102278,
                                                                       36152, 36182, 68282,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103391, 0, 3,
                                                                       102278, 67532, 102299,
                                                                       36182, 36212, 68327,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103454, 0, 3,
                                                                       102299, 67547, 102320,
                                                                       36212, 36242, 68372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103517, 0, 3,
                                                                       102320, 67562, 102341,
                                                                       36242, 36272, 68417,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103580, 0, 3,
                                                                       102341, 67577, 102362,
                                                                       36272, 36302, 68462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103643, 0, 3,
                                                                       102362, 67592, 102383,
                                                                       36302, 36332, 68507,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103706, 0, 3,
                                                                       102383, 67607, 102404,
                                                                       36332, 36362, 68552,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 103769, 0, 3,
                                                                       102404, 67622, 102425,
                                                                       36362, 36392, 68597,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 103832, 0, 3,
                                                                       102446, 67652, 102509,
                                                                       36452, 36512, 68642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 103958, 0, 3,
                                                                       102509, 67697, 102572,
                                                                       36512, 36572, 68732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 104084, 0, 3,
                                                                       102572, 67742, 102635,
                                                                       36572, 36632, 68822,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 104210, 0, 3,
                                                                       102635, 67787, 102698,
                                                                       36632, 36692, 68912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 104336, 0, 3,
                                                                       102698, 67832, 102761,
                                                                       36692, 36752, 69002,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 104462, 0, 3,
                                                                       102761, 67877, 102824,
                                                                       36752, 36812, 69092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 104588, 0, 3,
                                                                       102824, 67922, 102887,
                                                                       36812, 36872, 69182,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 104714, 0, 3,
                                                                       102887, 67967, 102950,
                                                                       36872, 36932, 69272,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 104840, 0, 3,
                                                                       102950, 68012, 103013,
                                                                       36932, 36992, 69362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 104966, 0, 3,
                                                                       103013, 68057, 103076,
                                                                       36992, 37052, 69452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 105092, 0, 3,
                                                                       103139, 68147, 103202,
                                                                       37172, 37232, 69542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 105218, 0, 3,
                                                                       103202, 68192, 103265,
                                                                       37232, 37292, 69632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 105344, 0, 3,
                                                                       103265, 68237, 103328,
                                                                       37292, 37352, 69722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 105470, 0, 3,
                                                                       103328, 68282, 103391,
                                                                       37352, 37412, 69812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 105596, 0, 3,
                                                                       103391, 68327, 103454,
                                                                       37412, 37472, 69902,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 105722, 0, 3,
                                                                       103454, 68372, 103517,
                                                                       37472, 37532, 69992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 105848, 0, 3,
                                                                       103517, 68417, 103580,
                                                                       37532, 37592, 70082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 105974, 0, 3,
                                                                       103580, 68462, 103643,
                                                                       37592, 37652, 70172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 106100, 0, 3,
                                                                       103643, 68507, 103706,
                                                                       37652, 37712, 70262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 106226, 0, 3,
                                                                       103706, 68552, 103769,
                                                                       37712, 37772, 70352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 106352, 0, 3,
                                                                       103832, 68642, 103958,
                                                                       37892, 37992, 70442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 106562, 0, 3,
                                                                       103958, 68732, 104084,
                                                                       37992, 38092, 70592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 106772, 0, 3,
                                                                       104084, 68822, 104210,
                                                                       38092, 38192, 70742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 106982, 0, 3,
                                                                       104210, 68912, 104336,
                                                                       38192, 38292, 70892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 107192, 0, 3,
                                                                       104336, 69002, 104462,
                                                                       38292, 38392, 71042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 107402, 0, 3,
                                                                       104462, 69092, 104588,
                                                                       38392, 38492, 71192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 107612, 0, 3,
                                                                       104588, 69182, 104714,
                                                                       38492, 38592, 71342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 107822, 0, 3,
                                                                       104714, 69272, 104840,
                                                                       38592, 38692, 71492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 108032, 0, 3,
                                                                       104840, 69362, 104966,
                                                                       38692, 38792, 71642,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 108242, 0, 3,
                                                                       105092, 69542, 105218,
                                                                       38992, 39092, 71792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 108452, 0, 3,
                                                                       105218, 69632, 105344,
                                                                       39092, 39192, 71942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 108662, 0, 3,
                                                                       105344, 69722, 105470,
                                                                       39192, 39292, 72092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 108872, 0, 3,
                                                                       105470, 69812, 105596,
                                                                       39292, 39392, 72242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109082, 0, 3,
                                                                       105596, 69902, 105722,
                                                                       39392, 39492, 72392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109292, 0, 3,
                                                                       105722, 69992, 105848,
                                                                       39492, 39592, 72542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109502, 0, 3,
                                                                       105848, 70082, 105974,
                                                                       39592, 39692, 72692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109712, 0, 3,
                                                                       105974, 70172, 106100,
                                                                       39692, 39792, 72842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109922, 0, 3,
                                                                       106100, 70262, 106226,
                                                                       39792, 39892, 72992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 110132, 0, 3,
                                                                       106352, 70442, 106562,
                                                                       40092, 40242, 73142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 110447, 0, 3,
                                                                       106562, 70592, 106772,
                                                                       40242, 40392, 73367,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 110762, 0, 3,
                                                                       106772, 70742, 106982,
                                                                       40392, 40542, 73592,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 111077, 0, 3,
                                                                       106982, 70892, 107192,
                                                                       40542, 40692, 73817,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 111392, 0, 3,
                                                                       107192, 71042, 107402,
                                                                       40692, 40842, 74042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 111707, 0, 3,
                                                                       107402, 71192, 107612,
                                                                       40842, 40992, 74267,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112022, 0, 3,
                                                                       107612, 71342, 107822,
                                                                       40992, 41142, 74492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112337, 0, 3,
                                                                       107822, 71492, 108032,
                                                                       41142, 41292, 74717,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112652, 0, 3,
                                                                       108242, 71792, 108452,
                                                                       41592, 41742, 74942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112967, 0, 3,
                                                                       108452, 71942, 108662,
                                                                       41742, 41892, 75167,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 113282, 0, 3,
                                                                       108662, 72092, 108872,
                                                                       41892, 42042, 75392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 113597, 0, 3,
                                                                       108872, 72242, 109082,
                                                                       42042, 42192, 75617,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 113912, 0, 3,
                                                                       109082, 72392, 109292,
                                                                       42192, 42342, 75842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 114227, 0, 3,
                                                                       109292, 72542, 109502,
                                                                       42342, 42492, 76067,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 114542, 0, 3,
                                                                       109502, 72692, 109712,
                                                                       42492, 42642, 76292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 114857, 0, 3,
                                                                       109712, 72842, 109922,
                                                                       42642, 42792, 76517,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 115172, 0, 3,
                                                                       110132, 73142, 110447,
                                                                       43092, 43302, 76742,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 115613, 0, 3,
                                                                       110447, 73367, 110762,
                                                                       43302, 43512, 77057,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 116054, 0, 3,
                                                                       110762, 73592, 111077,
                                                                       43512, 43722, 77372,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 116495, 0, 3,
                                                                       111077, 73817, 111392,
                                                                       43722, 43932, 77687,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 116936, 0, 3,
                                                                       111392, 74042, 111707,
                                                                       43932, 44142, 78002,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 117377, 0, 3,
                                                                       111707, 74267, 112022,
                                                                       44142, 44352, 78317,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 117818, 0, 3,
                                                                       112022, 74492, 112337,
                                                                       44352, 44562, 78632,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 118259, 0, 3,
                                                                       112652, 74942, 112967,
                                                                       44982, 45192, 78947,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 118700, 0, 3,
                                                                       112967, 75167, 113282,
                                                                       45192, 45402, 79262,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 119141, 0, 3,
                                                                       113282, 75392, 113597,
                                                                       45402, 45612, 79577,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 119582, 0, 3,
                                                                       113597, 75617, 113912,
                                                                       45612, 45822, 79892,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 120023, 0, 3,
                                                                       113912, 75842, 114227,
                                                                       45822, 46032, 80207,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 120464, 0, 3,
                                                                       114227, 76067, 114542,
                                                                       46032, 46242, 80522,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 120905, 0, 3,
                                                                       114542, 76292, 114857,
                                                                       46242, 46452, 80837,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 121346, 0, 3,
                                                                       115172, 76742, 115613,
                                                                       46872, 47152, 81152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 121934, 0, 3,
                                                                       115613, 77057, 116054,
                                                                       47152, 47432, 81572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 122522, 0, 3,
                                                                       116054, 77372, 116495,
                                                                       47432, 47712, 81992,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 123110, 0, 3,
                                                                       116495, 77687, 116936,
                                                                       47712, 47992, 82412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 123698, 0, 3,
                                                                       116936, 78002, 117377,
                                                                       47992, 48272, 82832,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 124286, 0, 3,
                                                                       117377, 78317, 117818,
                                                                       48272, 48552, 83252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 124874, 0, 3,
                                                                       118259, 78947, 118700,
                                                                       49112, 49392, 83672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 125462, 0, 3,
                                                                       118700, 79262, 119141,
                                                                       49392, 49672, 84092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 126050, 0, 3,
                                                                       119141, 79577, 119582,
                                                                       49672, 49952, 84512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 126638, 0, 3,
                                                                       119582, 79892, 120023,
                                                                       49952, 50232, 84932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 127226, 0, 3,
                                                                       120023, 80207, 120464,
                                                                       50232, 50512, 85352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 127814, 0, 3,
                                                                       120464, 80522, 120905,
                                                                       50512, 50792, 85772,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 128402, 0, 3,
                                                                       121346, 81152, 121934,
                                                                       51352, 51712, 86192,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 129158, 0, 3,
                                                                       121934, 81572, 122522,
                                                                       51712, 52072, 86732,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 129914, 0, 3,
                                                                       122522, 81992, 123110,
                                                                       52072, 52432, 87272,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 130670, 0, 3,
                                                                       123110, 82412, 123698,
                                                                       52432, 52792, 87812,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 131426, 0, 3,
                                                                       123698, 82832, 124286,
                                                                       52792, 53152, 88352,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 132182, 0, 3,
                                                                       124874, 83672, 125462,
                                                                       53872, 54232, 88892,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 132938, 0, 3,
                                                                       125462, 84092, 126050,
                                                                       54232, 54592, 89432,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 133694, 0, 3,
                                                                       126050, 84512, 126638,
                                                                       54592, 54952, 89972,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 134450, 0, 3,
                                                                       126638, 84932, 127226,
                                                                       54952, 55312, 90512,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 135206, 0, 3,
                                                                       127226, 85352, 127814,
                                                                       55312, 55672, 91052,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 135962, 0, 3,
                                                                       128402, 86192, 129158,
                                                                       56392, 56842, 91592,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 136907, 0, 3,
                                                                       129158, 86732, 129914,
                                                                       56842, 57292, 92267,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 137852, 0, 3,
                                                                       129914, 87272, 130670,
                                                                       57292, 57742, 92942,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 138797, 0, 3,
                                                                       130670, 87812, 131426,
                                                                       57742, 58192, 93617,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 139742, 0, 3,
                                                                       132182, 88892, 132938,
                                                                       59092, 59542, 94292,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 140687, 0, 3,
                                                                       132938, 89432, 133694,
                                                                       59542, 59992, 94967,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 141632, 0, 3,
                                                                       133694, 89972, 134450,
                                                                       59992, 60442, 95642,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 142577, 0, 3,
                                                                       134450, 90512, 135206,
                                                                       60442, 60892, 96317,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 143522, 0, 3,
                                                                       135962, 91592, 136907,
                                                                       61792, 62342, 96992,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 144677, 0, 3,
                                                                       136907, 92267, 137852,
                                                                       62342, 62892, 97817,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 145832, 0, 3,
                                                                       137852, 92942, 138797,
                                                                       62892, 63442, 98642,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 146987, 0, 3,
                                                                       139742, 94292, 140687,
                                                                       64542, 65092, 99467,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 148142, 0, 3,
                                                                       140687, 94967, 141632,
                                                                       65092, 65642, 100292,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 149297, 0, 3,
                                                                       141632, 95642, 142577,
                                                                       65642, 66192, 101117,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150452, 3, 67292,
                                                                       67307, 101984, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150480, 3, 67307,
                                                                       67322, 102005, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150508, 3, 67322,
                                                                       67337, 102026, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150536, 3, 67337,
                                                                       67352, 102047, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150564, 3, 67352,
                                                                       67367, 102068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150592, 3, 67367,
                                                                       67382, 102089, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150620, 3, 67382,
                                                                       67397, 102110, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150648, 3, 67397,
                                                                       67412, 102131, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150676, 3, 67412,
                                                                       67427, 102152, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150704, 3, 67427,
                                                                       67442, 102173, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150732, 3, 67472,
                                                                       67487, 102236, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150760, 3, 67487,
                                                                       67502, 102257, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150788, 3, 67502,
                                                                       67517, 102278, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150816, 3, 67517,
                                                                       67532, 102299, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150844, 3, 67532,
                                                                       67547, 102320, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150872, 3, 67547,
                                                                       67562, 102341, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150900, 3, 67562,
                                                                       67577, 102362, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150928, 3, 67577,
                                                                       67592, 102383, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150956, 3, 67592,
                                                                       67607, 102404, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150984, 3, 67607,
                                                                       67622, 102425, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151012, 0, 3,
                                                                       150452, 101984, 150480,
                                                                       67652, 67697, 102572,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151096, 0, 3,
                                                                       150480, 102005, 150508,
                                                                       67697, 67742, 102635,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151180, 0, 3,
                                                                       150508, 102026, 150536,
                                                                       67742, 67787, 102698,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151264, 0, 3,
                                                                       150536, 102047, 150564,
                                                                       67787, 67832, 102761,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151348, 0, 3,
                                                                       150564, 102068, 150592,
                                                                       67832, 67877, 102824,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151432, 0, 3,
                                                                       150592, 102089, 150620,
                                                                       67877, 67922, 102887,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151516, 0, 3,
                                                                       150620, 102110, 150648,
                                                                       67922, 67967, 102950,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151600, 0, 3,
                                                                       150648, 102131, 150676,
                                                                       67967, 68012, 103013,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151684, 0, 3,
                                                                       150676, 102152, 150704,
                                                                       68012, 68057, 103076,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151768, 0, 3,
                                                                       150732, 102236, 150760,
                                                                       68147, 68192, 103265,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151852, 0, 3,
                                                                       150760, 102257, 150788,
                                                                       68192, 68237, 103328,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151936, 0, 3,
                                                                       150788, 102278, 150816,
                                                                       68237, 68282, 103391,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 152020, 0, 3,
                                                                       150816, 102299, 150844,
                                                                       68282, 68327, 103454,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 152104, 0, 3,
                                                                       150844, 102320, 150872,
                                                                       68327, 68372, 103517,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 152188, 0, 3,
                                                                       150872, 102341, 150900,
                                                                       68372, 68417, 103580,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 152272, 0, 3,
                                                                       150900, 102362, 150928,
                                                                       68417, 68462, 103643,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 152356, 0, 3,
                                                                       150928, 102383, 150956,
                                                                       68462, 68507, 103706,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 152440, 0, 3,
                                                                       150956, 102404, 150984,
                                                                       68507, 68552, 103769,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152524, 0, 3,
                                                                       151012, 102572, 151096,
                                                                       68642, 68732, 104084,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152692, 0, 3,
                                                                       151096, 102635, 151180,
                                                                       68732, 68822, 104210,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152860, 0, 3,
                                                                       151180, 102698, 151264,
                                                                       68822, 68912, 104336,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153028, 0, 3,
                                                                       151264, 102761, 151348,
                                                                       68912, 69002, 104462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153196, 0, 3,
                                                                       151348, 102824, 151432,
                                                                       69002, 69092, 104588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153364, 0, 3,
                                                                       151432, 102887, 151516,
                                                                       69092, 69182, 104714,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153532, 0, 3,
                                                                       151516, 102950, 151600,
                                                                       69182, 69272, 104840,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153700, 0, 3,
                                                                       151600, 103013, 151684,
                                                                       69272, 69362, 104966,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153868, 0, 3,
                                                                       151768, 103265, 151852,
                                                                       69542, 69632, 105344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 154036, 0, 3,
                                                                       151852, 103328, 151936,
                                                                       69632, 69722, 105470,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 154204, 0, 3,
                                                                       151936, 103391, 152020,
                                                                       69722, 69812, 105596,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 154372, 0, 3,
                                                                       152020, 103454, 152104,
                                                                       69812, 69902, 105722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 154540, 0, 3,
                                                                       152104, 103517, 152188,
                                                                       69902, 69992, 105848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 154708, 0, 3,
                                                                       152188, 103580, 152272,
                                                                       69992, 70082, 105974,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 154876, 0, 3,
                                                                       152272, 103643, 152356,
                                                                       70082, 70172, 106100,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 155044, 0, 3,
                                                                       152356, 103706, 152440,
                                                                       70172, 70262, 106226,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155212, 0, 3,
                                                                       152524, 104084, 152692,
                                                                       70442, 70592, 106772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155492, 0, 3,
                                                                       152692, 104210, 152860,
                                                                       70592, 70742, 106982,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155772, 0, 3,
                                                                       152860, 104336, 153028,
                                                                       70742, 70892, 107192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 156052, 0, 3,
                                                                       153028, 104462, 153196,
                                                                       70892, 71042, 107402,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 156332, 0, 3,
                                                                       153196, 104588, 153364,
                                                                       71042, 71192, 107612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 156612, 0, 3,
                                                                       153364, 104714, 153532,
                                                                       71192, 71342, 107822,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 156892, 0, 3,
                                                                       153532, 104840, 153700,
                                                                       71342, 71492, 108032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 157172, 0, 3,
                                                                       153868, 105344, 154036,
                                                                       71792, 71942, 108662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 157452, 0, 3,
                                                                       154036, 105470, 154204,
                                                                       71942, 72092, 108872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 157732, 0, 3,
                                                                       154204, 105596, 154372,
                                                                       72092, 72242, 109082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 158012, 0, 3,
                                                                       154372, 105722, 154540,
                                                                       72242, 72392, 109292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 158292, 0, 3,
                                                                       154540, 105848, 154708,
                                                                       72392, 72542, 109502,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 158572, 0, 3,
                                                                       154708, 105974, 154876,
                                                                       72542, 72692, 109712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 158852, 0, 3,
                                                                       154876, 106100, 155044,
                                                                       72692, 72842, 109922,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 159132, 0, 3,
                                                                       155212, 106772, 155492,
                                                                       73142, 73367, 110762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 159552, 0, 3,
                                                                       155492, 106982, 155772,
                                                                       73367, 73592, 111077,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 159972, 0, 3,
                                                                       155772, 107192, 156052,
                                                                       73592, 73817, 111392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 160392, 0, 3,
                                                                       156052, 107402, 156332,
                                                                       73817, 74042, 111707,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 160812, 0, 3,
                                                                       156332, 107612, 156612,
                                                                       74042, 74267, 112022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 161232, 0, 3,
                                                                       156612, 107822, 156892,
                                                                       74267, 74492, 112337,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 161652, 0, 3,
                                                                       157172, 108662, 157452,
                                                                       74942, 75167, 113282,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 162072, 0, 3,
                                                                       157452, 108872, 157732,
                                                                       75167, 75392, 113597,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 162492, 0, 3,
                                                                       157732, 109082, 158012,
                                                                       75392, 75617, 113912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 162912, 0, 3,
                                                                       158012, 109292, 158292,
                                                                       75617, 75842, 114227,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 163332, 0, 3,
                                                                       158292, 109502, 158572,
                                                                       75842, 76067, 114542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 163752, 0, 3,
                                                                       158572, 109712, 158852,
                                                                       76067, 76292, 114857,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 164172, 0, 3,
                                                                       159132, 110762, 159552,
                                                                       76742, 77057, 116054,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 164760, 0, 3,
                                                                       159552, 111077, 159972,
                                                                       77057, 77372, 116495,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 165348, 0, 3,
                                                                       159972, 111392, 160392,
                                                                       77372, 77687, 116936,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 165936, 0, 3,
                                                                       160392, 111707, 160812,
                                                                       77687, 78002, 117377,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 166524, 0, 3,
                                                                       160812, 112022, 161232,
                                                                       78002, 78317, 117818,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 167112, 0, 3,
                                                                       161652, 113282, 162072,
                                                                       78947, 79262, 119141,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 167700, 0, 3,
                                                                       162072, 113597, 162492,
                                                                       79262, 79577, 119582,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 168288, 0, 3,
                                                                       162492, 113912, 162912,
                                                                       79577, 79892, 120023,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 168876, 0, 3,
                                                                       162912, 114227, 163332,
                                                                       79892, 80207, 120464,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 169464, 0, 3,
                                                                       163332, 114542, 163752,
                                                                       80207, 80522, 120905,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 170052, 0, 3,
                                                                       164172, 116054, 164760,
                                                                       81152, 81572, 122522,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 170836, 0, 3,
                                                                       164760, 116495, 165348,
                                                                       81572, 81992, 123110,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 171620, 0, 3,
                                                                       165348, 116936, 165936,
                                                                       81992, 82412, 123698,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 172404, 0, 3,
                                                                       165936, 117377, 166524,
                                                                       82412, 82832, 124286,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 173188, 0, 3,
                                                                       167112, 119141, 167700,
                                                                       83672, 84092, 126050,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 173972, 0, 3,
                                                                       167700, 119582, 168288,
                                                                       84092, 84512, 126638,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 174756, 0, 3,
                                                                       168288, 120023, 168876,
                                                                       84512, 84932, 127226,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 175540, 0, 3,
                                                                       168876, 120464, 169464,
                                                                       84932, 85352, 127814,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 176324, 0, 3,
                                                                       170052, 122522, 170836,
                                                                       86192, 86732, 129914,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 177332, 0, 3,
                                                                       170836, 123110, 171620,
                                                                       86732, 87272, 130670,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 178340, 0, 3,
                                                                       171620, 123698, 172404,
                                                                       87272, 87812, 131426,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 179348, 0, 3,
                                                                       173188, 126050, 173972,
                                                                       88892, 89432, 133694,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 180356, 0, 3,
                                                                       173972, 126638, 174756,
                                                                       89432, 89972, 134450,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 181364, 0, 3,
                                                                       174756, 127226, 175540,
                                                                       89972, 90512, 135206,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 182372, 0, 3,
                                                                       176324, 129914, 177332,
                                                                       91592, 92267, 137852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 183632, 0, 3,
                                                                       177332, 130670, 178340,
                                                                       92267, 92942, 138797,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 184892, 0, 3,
                                                                       179348, 133694, 180356,
                                                                       94292, 94967, 141632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 186152, 0, 3,
                                                                       180356, 134450, 181364,
                                                                       94967, 95642, 142577,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 187412, 0, 3,
                                                                       182372, 137852, 183632,
                                                                       96992, 97817, 145832,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 188952, 0, 3,
                                                                       184892, 141632, 186152,
                                                                       99467, 100292, 149297,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190492, 3, 101942,
                                                                       101963, 150452, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190528, 3, 101963,
                                                                       101984, 150480, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190564, 3, 101984,
                                                                       102005, 150508, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190600, 3, 102005,
                                                                       102026, 150536, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190636, 3, 102026,
                                                                       102047, 150564, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190672, 3, 102047,
                                                                       102068, 150592, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190708, 3, 102068,
                                                                       102089, 150620, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190744, 3, 102089,
                                                                       102110, 150648, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190780, 3, 102110,
                                                                       102131, 150676, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190816, 3, 102131,
                                                                       102152, 150704, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190852, 3, 102194,
                                                                       102215, 150732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190888, 3, 102215,
                                                                       102236, 150760, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190924, 3, 102236,
                                                                       102257, 150788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190960, 3, 102257,
                                                                       102278, 150816, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 190996, 3, 102278,
                                                                       102299, 150844, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191032, 3, 102299,
                                                                       102320, 150872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191068, 3, 102320,
                                                                       102341, 150900, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191104, 3, 102341,
                                                                       102362, 150928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191140, 3, 102362,
                                                                       102383, 150956, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 191176, 3, 102383,
                                                                       102404, 150984, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 191212, 0, 3,
                                                                       190492, 150452, 190528,
                                                                       102446, 102509, 151012,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 191320, 0, 3,
                                                                       190528, 150480, 190564,
                                                                       102509, 102572, 151096,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 191428, 0, 3,
                                                                       190564, 150508, 190600,
                                                                       102572, 102635, 151180,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 191536, 0, 3,
                                                                       190600, 150536, 190636,
                                                                       102635, 102698, 151264,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 191644, 0, 3,
                                                                       190636, 150564, 190672,
                                                                       102698, 102761, 151348,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 191752, 0, 3,
                                                                       190672, 150592, 190708,
                                                                       102761, 102824, 151432,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 191860, 0, 3,
                                                                       190708, 150620, 190744,
                                                                       102824, 102887, 151516,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 191968, 0, 3,
                                                                       190744, 150648, 190780,
                                                                       102887, 102950, 151600,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 192076, 0, 3,
                                                                       190780, 150676, 190816,
                                                                       102950, 103013, 151684,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 192184, 0, 3,
                                                                       190852, 150732, 190888,
                                                                       103139, 103202, 151768,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 192292, 0, 3,
                                                                       190888, 150760, 190924,
                                                                       103202, 103265, 151852,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 192400, 0, 3,
                                                                       190924, 150788, 190960,
                                                                       103265, 103328, 151936,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 192508, 0, 3,
                                                                       190960, 150816, 190996,
                                                                       103328, 103391, 152020,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 192616, 0, 3,
                                                                       190996, 150844, 191032,
                                                                       103391, 103454, 152104,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 192724, 0, 3,
                                                                       191032, 150872, 191068,
                                                                       103454, 103517, 152188,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 192832, 0, 3,
                                                                       191068, 150900, 191104,
                                                                       103517, 103580, 152272,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 192940, 0, 3,
                                                                       191104, 150928, 191140,
                                                                       103580, 103643, 152356,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 193048, 0, 3,
                                                                       191140, 150956, 191176,
                                                                       103643, 103706, 152440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 193156, 0, 3,
                                                                       191212, 151012, 191320,
                                                                       103832, 103958, 152524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 193372, 0, 3,
                                                                       191320, 151096, 191428,
                                                                       103958, 104084, 152692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 193588, 0, 3,
                                                                       191428, 151180, 191536,
                                                                       104084, 104210, 152860,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 193804, 0, 3,
                                                                       191536, 151264, 191644,
                                                                       104210, 104336, 153028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 194020, 0, 3,
                                                                       191644, 151348, 191752,
                                                                       104336, 104462, 153196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 194236, 0, 3,
                                                                       191752, 151432, 191860,
                                                                       104462, 104588, 153364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 194452, 0, 3,
                                                                       191860, 151516, 191968,
                                                                       104588, 104714, 153532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 194668, 0, 3,
                                                                       191968, 151600, 192076,
                                                                       104714, 104840, 153700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 194884, 0, 3,
                                                                       192184, 151768, 192292,
                                                                       105092, 105218, 153868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 195100, 0, 3,
                                                                       192292, 151852, 192400,
                                                                       105218, 105344, 154036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 195316, 0, 3,
                                                                       192400, 151936, 192508,
                                                                       105344, 105470, 154204,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 195532, 0, 3,
                                                                       192508, 152020, 192616,
                                                                       105470, 105596, 154372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 195748, 0, 3,
                                                                       192616, 152104, 192724,
                                                                       105596, 105722, 154540,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 195964, 0, 3,
                                                                       192724, 152188, 192832,
                                                                       105722, 105848, 154708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 196180, 0, 3,
                                                                       192832, 152272, 192940,
                                                                       105848, 105974, 154876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 196396, 0, 3,
                                                                       192940, 152356, 193048,
                                                                       105974, 106100, 155044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 196612, 0, 3,
                                                                       193156, 152524, 193372,
                                                                       106352, 106562, 155212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 196972, 0, 3,
                                                                       193372, 152692, 193588,
                                                                       106562, 106772, 155492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 197332, 0, 3,
                                                                       193588, 152860, 193804,
                                                                       106772, 106982, 155772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 197692, 0, 3,
                                                                       193804, 153028, 194020,
                                                                       106982, 107192, 156052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 198052, 0, 3,
                                                                       194020, 153196, 194236,
                                                                       107192, 107402, 156332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 198412, 0, 3,
                                                                       194236, 153364, 194452,
                                                                       107402, 107612, 156612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 198772, 0, 3,
                                                                       194452, 153532, 194668,
                                                                       107612, 107822, 156892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 199132, 0, 3,
                                                                       194884, 153868, 195100,
                                                                       108242, 108452, 157172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 199492, 0, 3,
                                                                       195100, 154036, 195316,
                                                                       108452, 108662, 157452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 199852, 0, 3,
                                                                       195316, 154204, 195532,
                                                                       108662, 108872, 157732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 200212, 0, 3,
                                                                       195532, 154372, 195748,
                                                                       108872, 109082, 158012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 200572, 0, 3,
                                                                       195748, 154540, 195964,
                                                                       109082, 109292, 158292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 200932, 0, 3,
                                                                       195964, 154708, 196180,
                                                                       109292, 109502, 158572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 201292, 0, 3,
                                                                       196180, 154876, 196396,
                                                                       109502, 109712, 158852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 201652, 0, 3,
                                                                       196612, 155212, 196972,
                                                                       110132, 110447, 159132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 202192, 0, 3,
                                                                       196972, 155492, 197332,
                                                                       110447, 110762, 159552,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 202732, 0, 3,
                                                                       197332, 155772, 197692,
                                                                       110762, 111077, 159972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 203272, 0, 3,
                                                                       197692, 156052, 198052,
                                                                       111077, 111392, 160392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 203812, 0, 3,
                                                                       198052, 156332, 198412,
                                                                       111392, 111707, 160812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 204352, 0, 3,
                                                                       198412, 156612, 198772,
                                                                       111707, 112022, 161232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 204892, 0, 3,
                                                                       199132, 157172, 199492,
                                                                       112652, 112967, 161652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 205432, 0, 3,
                                                                       199492, 157452, 199852,
                                                                       112967, 113282, 162072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 205972, 0, 3,
                                                                       199852, 157732, 200212,
                                                                       113282, 113597, 162492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 206512, 0, 3,
                                                                       200212, 158012, 200572,
                                                                       113597, 113912, 162912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 207052, 0, 3,
                                                                       200572, 158292, 200932,
                                                                       113912, 114227, 163332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 207592, 0, 3,
                                                                       200932, 158572, 201292,
                                                                       114227, 114542, 163752,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 208132, 0, 3,
                                                                       201652, 159132, 202192,
                                                                       115172, 115613, 164172,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 208888, 0, 3,
                                                                       202192, 159552, 202732,
                                                                       115613, 116054, 164760,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 209644, 0, 3,
                                                                       202732, 159972, 203272,
                                                                       116054, 116495, 165348,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 210400, 0, 3,
                                                                       203272, 160392, 203812,
                                                                       116495, 116936, 165936,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 211156, 0, 3,
                                                                       203812, 160812, 204352,
                                                                       116936, 117377, 166524,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 211912, 0, 3,
                                                                       204892, 161652, 205432,
                                                                       118259, 118700, 167112,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 212668, 0, 3,
                                                                       205432, 162072, 205972,
                                                                       118700, 119141, 167700,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 213424, 0, 3,
                                                                       205972, 162492, 206512,
                                                                       119141, 119582, 168288,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 214180, 0, 3,
                                                                       206512, 162912, 207052,
                                                                       119582, 120023, 168876,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 214936, 0, 3,
                                                                       207052, 163332, 207592,
                                                                       120023, 120464, 169464,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 215692, 0, 3,
                                                                       208132, 164172, 208888,
                                                                       121346, 121934, 170052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 216700, 0, 3,
                                                                       208888, 164760, 209644,
                                                                       121934, 122522, 170836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 217708, 0, 3,
                                                                       209644, 165348, 210400,
                                                                       122522, 123110, 171620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 218716, 0, 3,
                                                                       210400, 165936, 211156,
                                                                       123110, 123698, 172404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 219724, 0, 3,
                                                                       211912, 167112, 212668,
                                                                       124874, 125462, 173188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 220732, 0, 3,
                                                                       212668, 167700, 213424,
                                                                       125462, 126050, 173972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 221740, 0, 3,
                                                                       213424, 168288, 214180,
                                                                       126050, 126638, 174756,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 222748, 0, 3,
                                                                       214180, 168876, 214936,
                                                                       126638, 127226, 175540,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 223756, 0, 3,
                                                                       215692, 170052, 216700,
                                                                       128402, 129158, 176324,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 225052, 0, 3,
                                                                       216700, 170836, 217708,
                                                                       129158, 129914, 177332,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 226348, 0, 3,
                                                                       217708, 171620, 218716,
                                                                       129914, 130670, 178340,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 227644, 0, 3,
                                                                       219724, 173188, 220732,
                                                                       132182, 132938, 179348,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 228940, 0, 3,
                                                                       220732, 173972, 221740,
                                                                       132938, 133694, 180356,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 230236, 0, 3,
                                                                       221740, 174756, 222748,
                                                                       133694, 134450, 181364,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 231532, 0, 3,
                                                                       223756, 176324, 225052,
                                                                       135962, 136907, 182372,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 233152, 0, 3,
                                                                       225052, 177332, 226348,
                                                                       136907, 137852, 183632,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 234772, 0, 3,
                                                                       227644, 179348, 228940,
                                                                       139742, 140687, 184892,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 236392, 0, 3,
                                                                       228940, 180356, 230236,
                                                                       140687, 141632, 186152,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 238012, 0, 3,
                                                                       231532, 182372, 233152,
                                                                       143522, 144677, 187412,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 239992, 0, 3,
                                                                       234772, 184892, 236392,
                                                                       146987, 148142, 188952,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 241972, 208132, 756, ncols);

                    simdfunc::contract_primitives(buffer, 243043, 211912, 756, ncols);

                    simdfunc::contract_primitives(buffer, 244114, 215692, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 245542, 219724, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 246970, 223756, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 248806, 227644, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 250642, 231532, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 252937, 234772, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 255232, 238012, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 258037, 239992, 1980, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 242728, 241972, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 243799, 243043, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 245122, 244114, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 246550, 245542, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 248266, 246970, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 250102, 248806, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 252262, 250642, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 254557, 252937, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 257212, 255232, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 260017, 258037, 55, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 260842, 242728, 245122, 15, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 261787, 243799, 246550, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 262732, 245122, 248266, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 263992, 246550, 250102, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 265252, 248266, 252262, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 266872, 250102, 254557, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 268492, 252262, 257212, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 270517, 254557, 260017, 15, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 272542, 260842, 262732, 15, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 274432, 261787, 263992, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 276322, 262732, 265252, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 278842, 263992, 266872, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 281362, 265252, 268492, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 284602, 266872, 270517, 15, nmax);

        simdtrf::compute_hrr_fh(buffer, coordinates, 287842, 272542, 276322, 15, nmax);

        simdtrf::compute_hrr_fh(buffer, coordinates, 290992, 274432, 278842, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 294142, 276322, 281362, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 298342, 278842, 284602, 15, nmax);

        simdtrf::compute_hrr_gh(buffer, coordinates, 302542, 287842, 294142, 15, nmax);

        simdtrf::compute_hrr_gh(buffer, coordinates, 307267, 290992, 298342, 15, nmax);

        simdtrf::transform_h_inner(buffer, 311992, 307267, 15, 15, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 311992, 165, nmax);

        simdtrf::transform_h_inner(buffer, 311992, 302542, 15, 15, nmax);

        simdtrf::transform_g_outer(values + 1485 * nvalues + n * npairs, nvalues, buffer, 311992,
                                   165, nmax);
    }

    for (size_t m = 0; m < 2970; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
