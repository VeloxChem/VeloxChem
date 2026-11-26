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


#include "SimdThreeCenterElectronRepulsionRsRecHIH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"
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
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_hih_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hih_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 305375, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3146 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 305375, 187796, 18854, dimensions);

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

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4078, 0, 3, 2588,
                                                                       2633, 3308, 3363, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4144, 0, 3, 2633,
                                                                       2678, 3363, 3418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4210, 0, 3, 2678,
                                                                       2723, 3418, 3473, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4276, 0, 3, 2723,
                                                                       2768, 3473, 3528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4342, 0, 3, 2768,
                                                                       2813, 3528, 3583, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4408, 0, 3, 2813,
                                                                       2858, 3583, 3638, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4474, 0, 3, 2948,
                                                                       2993, 3693, 3748, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4540, 0, 3, 2993,
                                                                       3038, 3748, 3803, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4606, 0, 3, 3038,
                                                                       3083, 3803, 3858, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4672, 0, 3, 3083,
                                                                       3128, 3858, 3913, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4738, 0, 3, 3128,
                                                                       3173, 3913, 3968, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4804, 0, 3, 3173,
                                                                       3218, 3968, 4023, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 4870, 0, 3, 3308,
                                                                       3363, 4078, 4144, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 4948, 0, 3, 3363,
                                                                       3418, 4144, 4210, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 5026, 0, 3, 3418,
                                                                       3473, 4210, 4276, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 5104, 0, 3, 3473,
                                                                       3528, 4276, 4342, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 5182, 0, 3, 3528,
                                                                       3583, 4342, 4408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 5260, 0, 3, 3693,
                                                                       3748, 4474, 4540, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 5338, 0, 3, 3748,
                                                                       3803, 4540, 4606, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 5416, 0, 3, 3803,
                                                                       3858, 4606, 4672, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 5494, 0, 3, 3858,
                                                                       3913, 4672, 4738, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 5572, 0, 3, 3913,
                                                                       3968, 4738, 4804, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5650, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5653, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5656, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5659, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5662, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5665, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5668, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5671, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5674, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5677, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5680, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5683, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5686, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5689, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5692, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5695, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5698, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5701, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5704, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5707, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5710, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5713, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5716, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5719, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5722, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5725, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5728, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5731, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5734, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5737, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5740, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 5743, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5746, 3, 9, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5755, 3, 10, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5764, 3, 11, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5773, 3, 12, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5782, 3, 13, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5791, 3, 14, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5800, 3, 15, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5809, 3, 16, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5818, 3, 17, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5827, 3, 18, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5836, 3, 19, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5845, 3, 20, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5854, 3, 21, 82,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5863, 3, 26, 91,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5872, 3, 27, 94,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5881, 3, 28, 97,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5890, 3, 29, 100,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5899, 3, 30, 103,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5908, 3, 31, 106,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5917, 3, 32, 109,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5926, 3, 33, 112,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5935, 3, 34, 115,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5944, 3, 35, 118,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5953, 3, 36, 121,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5962, 3, 37, 124,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 5971, 3, 38, 127,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5980, 3, 40, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 5998, 3, 43, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6016, 3, 46, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6034, 3, 49, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6052, 3, 52, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6070, 3, 55, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6088, 3, 58, 166,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6106, 3, 61, 172,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6124, 3, 64, 178,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6142, 3, 67, 184,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6160, 3, 70, 190,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6178, 3, 73, 196,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6196, 3, 76, 202,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6214, 3, 79, 208,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6232, 3, 85, 214,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6250, 3, 88, 220,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6268, 3, 91, 226,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6286, 3, 94, 232,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6304, 3, 97, 238,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6322, 3, 100, 244,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6340, 3, 103, 250,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6358, 3, 106, 256,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6376, 3, 109, 262,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6394, 3, 112, 268,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6412, 3, 115, 274,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6430, 3, 118, 280,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6448, 3, 121, 286,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 6466, 3, 124, 292,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6484, 3, 130, 298,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6514, 3, 136, 308,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6544, 3, 142, 318,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6574, 3, 148, 328,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6604, 3, 154, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6634, 3, 160, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6664, 3, 166, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6694, 3, 172, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6724, 3, 178, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6754, 3, 184, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6784, 3, 190, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6814, 3, 196, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6844, 3, 202, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6874, 3, 214, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6904, 3, 220, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6934, 3, 226, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6964, 3, 232, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 6994, 3, 238, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7024, 3, 244, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7054, 3, 250, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7084, 3, 256, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7114, 3, 262, 508,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7144, 3, 268, 518,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7174, 3, 274, 528,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7204, 3, 280, 538,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 7234, 3, 286, 548,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7264, 3, 298, 558,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7309, 3, 308, 573,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7354, 3, 318, 588,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7399, 3, 328, 603,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7444, 3, 338, 618,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7489, 3, 348, 633,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7534, 3, 358, 648,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7579, 3, 368, 663,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7624, 3, 378, 678,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7669, 3, 388, 693,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7714, 3, 398, 708,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7759, 3, 408, 723,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7804, 3, 428, 738,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7849, 3, 438, 753,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7894, 3, 448, 768,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7939, 3, 458, 783,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 7984, 3, 468, 798,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8029, 3, 478, 813,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8074, 3, 488, 828,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8119, 3, 498, 843,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8164, 3, 508, 858,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8209, 3, 518, 873,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8254, 3, 528, 888,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8299, 3, 538, 903,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8344, 3, 558, 918,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8407, 3, 573, 939,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8470, 3, 588, 960,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8533, 3, 603, 981,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8596, 3, 618,
                                                                       1002, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8659, 3, 633,
                                                                       1023, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8722, 3, 648,
                                                                       1044, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8785, 3, 663,
                                                                       1065, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8848, 3, 678,
                                                                       1086, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8911, 3, 693,
                                                                       1107, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 8974, 3, 708,
                                                                       1128, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9037, 3, 738,
                                                                       1149, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9100, 3, 753,
                                                                       1170, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9163, 3, 768,
                                                                       1191, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9226, 3, 783,
                                                                       1212, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9289, 3, 798,
                                                                       1233, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9352, 3, 813,
                                                                       1254, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9415, 3, 828,
                                                                       1275, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9478, 3, 843,
                                                                       1296, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9541, 3, 858,
                                                                       1317, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9604, 3, 873,
                                                                       1338, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 9667, 3, 888,
                                                                       1359, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9730, 3, 918,
                                                                       1380, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9814, 3, 939,
                                                                       1408, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9898, 3, 960,
                                                                       1436, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 9982, 3, 981,
                                                                       1464, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10066, 3, 1002,
                                                                       1492, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10150, 3, 1023,
                                                                       1520, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10234, 3, 1044,
                                                                       1548, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10318, 3, 1065,
                                                                       1576, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10402, 3, 1086,
                                                                       1604, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10486, 3, 1107,
                                                                       1632, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10570, 3, 1149,
                                                                       1660, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10654, 3, 1170,
                                                                       1688, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10738, 3, 1191,
                                                                       1716, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10822, 3, 1212,
                                                                       1744, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10906, 3, 1233,
                                                                       1772, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 10990, 3, 1254,
                                                                       1800, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11074, 3, 1275,
                                                                       1828, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11158, 3, 1296,
                                                                       1856, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11242, 3, 1317,
                                                                       1884, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11326, 3, 1338,
                                                                       1912, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11410, 3, 1380,
                                                                       1940, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11518, 3, 1408,
                                                                       1976, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11626, 3, 1436,
                                                                       2012, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11734, 3, 1464,
                                                                       2048, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11842, 3, 1492,
                                                                       2084, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 11950, 3, 1520,
                                                                       2120, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12058, 3, 1548,
                                                                       2156, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12166, 3, 1576,
                                                                       2192, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12274, 3, 1604,
                                                                       2228, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12382, 3, 1660,
                                                                       2264, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12490, 3, 1688,
                                                                       2300, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12598, 3, 1716,
                                                                       2336, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12706, 3, 1744,
                                                                       2372, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12814, 3, 1772,
                                                                       2408, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 12922, 3, 1800,
                                                                       2444, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13030, 3, 1828,
                                                                       2480, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13138, 3, 1856,
                                                                       2516, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13246, 3, 1884,
                                                                       2552, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13354, 3, 1940,
                                                                       2588, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13489, 3, 1976,
                                                                       2633, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13624, 3, 2012,
                                                                       2678, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13759, 3, 2048,
                                                                       2723, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 13894, 3, 2084,
                                                                       2768, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14029, 3, 2120,
                                                                       2813, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14164, 3, 2156,
                                                                       2858, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14299, 3, 2192,
                                                                       2903, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14434, 3, 2264,
                                                                       2948, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14569, 3, 2300,
                                                                       2993, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14704, 3, 2336,
                                                                       3038, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14839, 3, 2372,
                                                                       3083, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 14974, 3, 2408,
                                                                       3128, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15109, 3, 2444,
                                                                       3173, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15244, 3, 2480,
                                                                       3218, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 15379, 3, 2516,
                                                                       3263, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 15514, 3, 2588,
                                                                       3308, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 15679, 3, 2633,
                                                                       3363, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 15844, 3, 2678,
                                                                       3418, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16009, 3, 2723,
                                                                       3473, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16174, 3, 2768,
                                                                       3528, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16339, 3, 2813,
                                                                       3583, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16504, 3, 2858,
                                                                       3638, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16669, 3, 2948,
                                                                       3693, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16834, 3, 2993,
                                                                       3748, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 16999, 3, 3038,
                                                                       3803, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17164, 3, 3083,
                                                                       3858, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17329, 3, 3128,
                                                                       3913, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17494, 3, 3173,
                                                                       3968, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 17659, 3, 3218,
                                                                       4023, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 17824, 3, 3308,
                                                                       4078, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 18022, 3, 3363,
                                                                       4144, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 18220, 3, 3418,
                                                                       4210, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 18418, 3, 3473,
                                                                       4276, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 18616, 3, 3528,
                                                                       4342, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 18814, 3, 3583,
                                                                       4408, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 19012, 3, 3693,
                                                                       4474, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 19210, 3, 3748,
                                                                       4540, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 19408, 3, 3803,
                                                                       4606, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 19606, 3, 3858,
                                                                       4672, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 19804, 3, 3913,
                                                                       4738, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 20002, 3, 3968,
                                                                       4804, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 20200, 3, 4078,
                                                                       4870, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 20434, 3, 4144,
                                                                       4948, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 20668, 3, 4210,
                                                                       5026, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 20902, 3, 4276,
                                                                       5104, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 21136, 3, 4342,
                                                                       5182, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 21370, 3, 4474,
                                                                       5260, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 21604, 3, 4540,
                                                                       5338, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 21838, 3, 4606,
                                                                       5416, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 22072, 3, 4672,
                                                                       5494, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 22306, 3, 4738,
                                                                       5572, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22540, 3, 7, 8,
                                                                       5656, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22546, 3, 8, 9,
                                                                       5659, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22552, 3, 9, 10,
                                                                       5662, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22558, 3, 10, 11,
                                                                       5665, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22564, 3, 11, 12,
                                                                       5668, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22570, 3, 12, 13,
                                                                       5671, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22576, 3, 13, 14,
                                                                       5674, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22582, 3, 14, 15,
                                                                       5677, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22588, 3, 15, 16,
                                                                       5680, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22594, 3, 16, 17,
                                                                       5683, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22600, 3, 17, 18,
                                                                       5686, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22606, 3, 18, 19,
                                                                       5689, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22612, 3, 19, 20,
                                                                       5692, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22618, 3, 20, 21,
                                                                       5695, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22624, 3, 24, 25,
                                                                       5704, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22630, 3, 25, 26,
                                                                       5707, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22636, 3, 26, 27,
                                                                       5710, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22642, 3, 27, 28,
                                                                       5713, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22648, 3, 28, 29,
                                                                       5716, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22654, 3, 29, 30,
                                                                       5719, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22660, 3, 30, 31,
                                                                       5722, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22666, 3, 31, 32,
                                                                       5725, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22672, 3, 32, 33,
                                                                       5728, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22678, 3, 33, 34,
                                                                       5731, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22684, 3, 34, 35,
                                                                       5734, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22690, 3, 35, 36,
                                                                       5737, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22696, 3, 36, 37,
                                                                       5740, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 22702, 3, 37, 38,
                                                                       5743, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22708, 0, 3,
                                                                       22540, 5656, 22546, 5746,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22726, 0, 3,
                                                                       22546, 5659, 22552, 5755,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22744, 0, 3,
                                                                       22552, 5662, 22558, 5764,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22762, 0, 3,
                                                                       22558, 5665, 22564, 5773,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22780, 0, 3,
                                                                       22564, 5668, 22570, 5782,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22798, 0, 3,
                                                                       22570, 5671, 22576, 5791,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22816, 0, 3,
                                                                       22576, 5674, 22582, 5800,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22834, 0, 3,
                                                                       22582, 5677, 22588, 5809,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22852, 0, 3,
                                                                       22588, 5680, 22594, 5818,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22870, 0, 3,
                                                                       22594, 5683, 22600, 5827,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22888, 0, 3,
                                                                       22600, 5686, 22606, 5836,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22906, 0, 3,
                                                                       22606, 5689, 22612, 5845,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22924, 0, 3,
                                                                       22612, 5692, 22618, 5854,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22942, 0, 3,
                                                                       22624, 5704, 22630, 5863,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22960, 0, 3,
                                                                       22630, 5707, 22636, 5872,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22978, 0, 3,
                                                                       22636, 5710, 22642, 5881,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 22996, 0, 3,
                                                                       22642, 5713, 22648, 5890,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23014, 0, 3,
                                                                       22648, 5716, 22654, 5899,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23032, 0, 3,
                                                                       22654, 5719, 22660, 5908,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23050, 0, 3,
                                                                       22660, 5722, 22666, 5917,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23068, 0, 3,
                                                                       22666, 5725, 22672, 5926,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23086, 0, 3,
                                                                       22672, 5728, 22678, 5935,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23104, 0, 3,
                                                                       22678, 5731, 22684, 5944,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23122, 0, 3,
                                                                       22684, 5734, 22690, 5953,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23140, 0, 3,
                                                                       22690, 5737, 22696, 5962,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 23158, 0, 3,
                                                                       22696, 5740, 22702, 5971,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23176, 0, 3,
                                                                       22708, 5746, 22726, 130,
                                                                       136, 6016, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23212, 0, 3,
                                                                       22726, 5755, 22744, 136,
                                                                       142, 6034, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23248, 0, 3,
                                                                       22744, 5764, 22762, 142,
                                                                       148, 6052, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23284, 0, 3,
                                                                       22762, 5773, 22780, 148,
                                                                       154, 6070, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23320, 0, 3,
                                                                       22780, 5782, 22798, 154,
                                                                       160, 6088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23356, 0, 3,
                                                                       22798, 5791, 22816, 160,
                                                                       166, 6106, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23392, 0, 3,
                                                                       22816, 5800, 22834, 166,
                                                                       172, 6124, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23428, 0, 3,
                                                                       22834, 5809, 22852, 172,
                                                                       178, 6142, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23464, 0, 3,
                                                                       22852, 5818, 22870, 178,
                                                                       184, 6160, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23500, 0, 3,
                                                                       22870, 5827, 22888, 184,
                                                                       190, 6178, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23536, 0, 3,
                                                                       22888, 5836, 22906, 190,
                                                                       196, 6196, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23572, 0, 3,
                                                                       22906, 5845, 22924, 196,
                                                                       202, 6214, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23608, 0, 3,
                                                                       22942, 5863, 22960, 214,
                                                                       220, 6268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23644, 0, 3,
                                                                       22960, 5872, 22978, 220,
                                                                       226, 6286, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23680, 0, 3,
                                                                       22978, 5881, 22996, 226,
                                                                       232, 6304, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23716, 0, 3,
                                                                       22996, 5890, 23014, 232,
                                                                       238, 6322, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23752, 0, 3,
                                                                       23014, 5899, 23032, 238,
                                                                       244, 6340, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23788, 0, 3,
                                                                       23032, 5908, 23050, 244,
                                                                       250, 6358, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23824, 0, 3,
                                                                       23050, 5917, 23068, 250,
                                                                       256, 6376, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23860, 0, 3,
                                                                       23068, 5926, 23086, 256,
                                                                       262, 6394, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23896, 0, 3,
                                                                       23086, 5935, 23104, 262,
                                                                       268, 6412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23932, 0, 3,
                                                                       23104, 5944, 23122, 268,
                                                                       274, 6430, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 23968, 0, 3,
                                                                       23122, 5953, 23140, 274,
                                                                       280, 6448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 24004, 0, 3,
                                                                       23140, 5962, 23158, 280,
                                                                       286, 6466, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24040, 0, 3,
                                                                       23176, 6016, 23212, 298,
                                                                       308, 6544, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24100, 0, 3,
                                                                       23212, 6034, 23248, 308,
                                                                       318, 6574, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24160, 0, 3,
                                                                       23248, 6052, 23284, 318,
                                                                       328, 6604, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24220, 0, 3,
                                                                       23284, 6070, 23320, 328,
                                                                       338, 6634, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24280, 0, 3,
                                                                       23320, 6088, 23356, 338,
                                                                       348, 6664, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24340, 0, 3,
                                                                       23356, 6106, 23392, 348,
                                                                       358, 6694, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24400, 0, 3,
                                                                       23392, 6124, 23428, 358,
                                                                       368, 6724, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24460, 0, 3,
                                                                       23428, 6142, 23464, 368,
                                                                       378, 6754, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24520, 0, 3,
                                                                       23464, 6160, 23500, 378,
                                                                       388, 6784, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24580, 0, 3,
                                                                       23500, 6178, 23536, 388,
                                                                       398, 6814, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24640, 0, 3,
                                                                       23536, 6196, 23572, 398,
                                                                       408, 6844, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24700, 0, 3,
                                                                       23608, 6268, 23644, 428,
                                                                       438, 6934, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24760, 0, 3,
                                                                       23644, 6286, 23680, 438,
                                                                       448, 6964, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24820, 0, 3,
                                                                       23680, 6304, 23716, 448,
                                                                       458, 6994, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24880, 0, 3,
                                                                       23716, 6322, 23752, 458,
                                                                       468, 7024, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 24940, 0, 3,
                                                                       23752, 6340, 23788, 468,
                                                                       478, 7054, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25000, 0, 3,
                                                                       23788, 6358, 23824, 478,
                                                                       488, 7084, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25060, 0, 3,
                                                                       23824, 6376, 23860, 488,
                                                                       498, 7114, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25120, 0, 3,
                                                                       23860, 6394, 23896, 498,
                                                                       508, 7144, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25180, 0, 3,
                                                                       23896, 6412, 23932, 508,
                                                                       518, 7174, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25240, 0, 3,
                                                                       23932, 6430, 23968, 518,
                                                                       528, 7204, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 25300, 0, 3,
                                                                       23968, 6448, 24004, 528,
                                                                       538, 7234, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25360, 0, 3,
                                                                       24040, 6544, 24100, 558,
                                                                       573, 7354, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25450, 0, 3,
                                                                       24100, 6574, 24160, 573,
                                                                       588, 7399, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25540, 0, 3,
                                                                       24160, 6604, 24220, 588,
                                                                       603, 7444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25630, 0, 3,
                                                                       24220, 6634, 24280, 603,
                                                                       618, 7489, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25720, 0, 3,
                                                                       24280, 6664, 24340, 618,
                                                                       633, 7534, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25810, 0, 3,
                                                                       24340, 6694, 24400, 633,
                                                                       648, 7579, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25900, 0, 3,
                                                                       24400, 6724, 24460, 648,
                                                                       663, 7624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 25990, 0, 3,
                                                                       24460, 6754, 24520, 663,
                                                                       678, 7669, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26080, 0, 3,
                                                                       24520, 6784, 24580, 678,
                                                                       693, 7714, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26170, 0, 3,
                                                                       24580, 6814, 24640, 693,
                                                                       708, 7759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26260, 0, 3,
                                                                       24700, 6934, 24760, 738,
                                                                       753, 7894, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26350, 0, 3,
                                                                       24760, 6964, 24820, 753,
                                                                       768, 7939, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26440, 0, 3,
                                                                       24820, 6994, 24880, 768,
                                                                       783, 7984, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26530, 0, 3,
                                                                       24880, 7024, 24940, 783,
                                                                       798, 8029, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26620, 0, 3,
                                                                       24940, 7054, 25000, 798,
                                                                       813, 8074, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26710, 0, 3,
                                                                       25000, 7084, 25060, 813,
                                                                       828, 8119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26800, 0, 3,
                                                                       25060, 7114, 25120, 828,
                                                                       843, 8164, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26890, 0, 3,
                                                                       25120, 7144, 25180, 843,
                                                                       858, 8209, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 26980, 0, 3,
                                                                       25180, 7174, 25240, 858,
                                                                       873, 8254, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 27070, 0, 3,
                                                                       25240, 7204, 25300, 873,
                                                                       888, 8299, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27160, 0, 3,
                                                                       25360, 7354, 25450, 918,
                                                                       939, 8470, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27286, 0, 3,
                                                                       25450, 7399, 25540, 939,
                                                                       960, 8533, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27412, 0, 3,
                                                                       25540, 7444, 25630, 960,
                                                                       981, 8596, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27538, 0, 3,
                                                                       25630, 7489, 25720, 981,
                                                                       1002, 8659, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27664, 0, 3,
                                                                       25720, 7534, 25810, 1002,
                                                                       1023, 8722, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27790, 0, 3,
                                                                       25810, 7579, 25900, 1023,
                                                                       1044, 8785, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 27916, 0, 3,
                                                                       25900, 7624, 25990, 1044,
                                                                       1065, 8848, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28042, 0, 3,
                                                                       25990, 7669, 26080, 1065,
                                                                       1086, 8911, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28168, 0, 3,
                                                                       26080, 7714, 26170, 1086,
                                                                       1107, 8974, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28294, 0, 3,
                                                                       26260, 7894, 26350, 1149,
                                                                       1170, 9163, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28420, 0, 3,
                                                                       26350, 7939, 26440, 1170,
                                                                       1191, 9226, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28546, 0, 3,
                                                                       26440, 7984, 26530, 1191,
                                                                       1212, 9289, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28672, 0, 3,
                                                                       26530, 8029, 26620, 1212,
                                                                       1233, 9352, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28798, 0, 3,
                                                                       26620, 8074, 26710, 1233,
                                                                       1254, 9415, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 28924, 0, 3,
                                                                       26710, 8119, 26800, 1254,
                                                                       1275, 9478, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29050, 0, 3,
                                                                       26800, 8164, 26890, 1275,
                                                                       1296, 9541, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29176, 0, 3,
                                                                       26890, 8209, 26980, 1296,
                                                                       1317, 9604, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 29302, 0, 3,
                                                                       26980, 8254, 27070, 1317,
                                                                       1338, 9667, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 29428, 0, 3,
                                                                       27160, 8470, 27286, 1380,
                                                                       1408, 9898, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 29596, 0, 3,
                                                                       27286, 8533, 27412, 1408,
                                                                       1436, 9982, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 29764, 0, 3,
                                                                       27412, 8596, 27538, 1436,
                                                                       1464, 10066, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 29932, 0, 3,
                                                                       27538, 8659, 27664, 1464,
                                                                       1492, 10150, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30100, 0, 3,
                                                                       27664, 8722, 27790, 1492,
                                                                       1520, 10234, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30268, 0, 3,
                                                                       27790, 8785, 27916, 1520,
                                                                       1548, 10318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30436, 0, 3,
                                                                       27916, 8848, 28042, 1548,
                                                                       1576, 10402, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30604, 0, 3,
                                                                       28042, 8911, 28168, 1576,
                                                                       1604, 10486, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30772, 0, 3,
                                                                       28294, 9163, 28420, 1660,
                                                                       1688, 10738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 30940, 0, 3,
                                                                       28420, 9226, 28546, 1688,
                                                                       1716, 10822, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31108, 0, 3,
                                                                       28546, 9289, 28672, 1716,
                                                                       1744, 10906, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31276, 0, 3,
                                                                       28672, 9352, 28798, 1744,
                                                                       1772, 10990, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31444, 0, 3,
                                                                       28798, 9415, 28924, 1772,
                                                                       1800, 11074, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31612, 0, 3,
                                                                       28924, 9478, 29050, 1800,
                                                                       1828, 11158, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31780, 0, 3,
                                                                       29050, 9541, 29176, 1828,
                                                                       1856, 11242, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 31948, 0, 3,
                                                                       29176, 9604, 29302, 1856,
                                                                       1884, 11326, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 32116, 0, 3,
                                                                       29428, 9898, 29596, 1940,
                                                                       1976, 11626, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 32332, 0, 3,
                                                                       29596, 9982, 29764, 1976,
                                                                       2012, 11734, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 32548, 0, 3,
                                                                       29764, 10066, 29932, 2012,
                                                                       2048, 11842, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 32764, 0, 3,
                                                                       29932, 10150, 30100, 2048,
                                                                       2084, 11950, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 32980, 0, 3,
                                                                       30100, 10234, 30268, 2084,
                                                                       2120, 12058, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 33196, 0, 3,
                                                                       30268, 10318, 30436, 2120,
                                                                       2156, 12166, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 33412, 0, 3,
                                                                       30436, 10402, 30604, 2156,
                                                                       2192, 12274, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 33628, 0, 3,
                                                                       30772, 10738, 30940, 2264,
                                                                       2300, 12598, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 33844, 0, 3,
                                                                       30940, 10822, 31108, 2300,
                                                                       2336, 12706, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 34060, 0, 3,
                                                                       31108, 10906, 31276, 2336,
                                                                       2372, 12814, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 34276, 0, 3,
                                                                       31276, 10990, 31444, 2372,
                                                                       2408, 12922, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 34492, 0, 3,
                                                                       31444, 11074, 31612, 2408,
                                                                       2444, 13030, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 34708, 0, 3,
                                                                       31612, 11158, 31780, 2444,
                                                                       2480, 13138, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 34924, 0, 3,
                                                                       31780, 11242, 31948, 2480,
                                                                       2516, 13246, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 35140, 0, 3,
                                                                       32116, 11626, 32332, 2588,
                                                                       2633, 13624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 35410, 0, 3,
                                                                       32332, 11734, 32548, 2633,
                                                                       2678, 13759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 35680, 0, 3,
                                                                       32548, 11842, 32764, 2678,
                                                                       2723, 13894, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 35950, 0, 3,
                                                                       32764, 11950, 32980, 2723,
                                                                       2768, 14029, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 36220, 0, 3,
                                                                       32980, 12058, 33196, 2768,
                                                                       2813, 14164, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 36490, 0, 3,
                                                                       33196, 12166, 33412, 2813,
                                                                       2858, 14299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 36760, 0, 3,
                                                                       33628, 12598, 33844, 2948,
                                                                       2993, 14704, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 37030, 0, 3,
                                                                       33844, 12706, 34060, 2993,
                                                                       3038, 14839, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 37300, 0, 3,
                                                                       34060, 12814, 34276, 3038,
                                                                       3083, 14974, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 37570, 0, 3,
                                                                       34276, 12922, 34492, 3083,
                                                                       3128, 15109, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 37840, 0, 3,
                                                                       34492, 13030, 34708, 3128,
                                                                       3173, 15244, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 38110, 0, 3,
                                                                       34708, 13138, 34924, 3173,
                                                                       3218, 15379, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 38380, 0, 3,
                                                                       35140, 13624, 35410, 3308,
                                                                       3363, 15844, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 38710, 0, 3,
                                                                       35410, 13759, 35680, 3363,
                                                                       3418, 16009, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 39040, 0, 3,
                                                                       35680, 13894, 35950, 3418,
                                                                       3473, 16174, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 39370, 0, 3,
                                                                       35950, 14029, 36220, 3473,
                                                                       3528, 16339, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 39700, 0, 3,
                                                                       36220, 14164, 36490, 3528,
                                                                       3583, 16504, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 40030, 0, 3,
                                                                       36760, 14704, 37030, 3693,
                                                                       3748, 16999, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 40360, 0, 3,
                                                                       37030, 14839, 37300, 3748,
                                                                       3803, 17164, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 40690, 0, 3,
                                                                       37300, 14974, 37570, 3803,
                                                                       3858, 17329, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 41020, 0, 3,
                                                                       37570, 15109, 37840, 3858,
                                                                       3913, 17494, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 41350, 0, 3,
                                                                       37840, 15244, 38110, 3913,
                                                                       3968, 17659, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 41680, 0, 3,
                                                                       38380, 15844, 38710, 4078,
                                                                       4144, 18220, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 42076, 0, 3,
                                                                       38710, 16009, 39040, 4144,
                                                                       4210, 18418, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 42472, 0, 3,
                                                                       39040, 16174, 39370, 4210,
                                                                       4276, 18616, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 42868, 0, 3,
                                                                       39370, 16339, 39700, 4276,
                                                                       4342, 18814, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 43264, 0, 3,
                                                                       40030, 16999, 40360, 4474,
                                                                       4540, 19408, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 43660, 0, 3,
                                                                       40360, 17164, 40690, 4540,
                                                                       4606, 19606, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 44056, 0, 3,
                                                                       40690, 17329, 41020, 4606,
                                                                       4672, 19804, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 44452, 0, 3,
                                                                       41020, 17494, 41350, 4672,
                                                                       4738, 20002, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 44848, 0, 3,
                                                                       41680, 18220, 42076, 4870,
                                                                       4948, 20668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 45316, 0, 3,
                                                                       42076, 18418, 42472, 4948,
                                                                       5026, 20902, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 45784, 0, 3,
                                                                       42472, 18616, 42868, 5026,
                                                                       5104, 21136, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 46252, 0, 3,
                                                                       43264, 19408, 43660, 5260,
                                                                       5338, 21838, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 46720, 0, 3,
                                                                       43660, 19606, 44056, 5338,
                                                                       5416, 22072, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 47188, 0, 3,
                                                                       44056, 19804, 44452, 5416,
                                                                       5494, 22306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47656, 3, 5650,
                                                                       5653, 22540, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47666, 3, 5653,
                                                                       5656, 22546, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47676, 3, 5656,
                                                                       5659, 22552, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47686, 3, 5659,
                                                                       5662, 22558, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47696, 3, 5662,
                                                                       5665, 22564, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47706, 3, 5665,
                                                                       5668, 22570, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47716, 3, 5668,
                                                                       5671, 22576, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47726, 3, 5671,
                                                                       5674, 22582, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47736, 3, 5674,
                                                                       5677, 22588, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47746, 3, 5677,
                                                                       5680, 22594, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47756, 3, 5680,
                                                                       5683, 22600, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47766, 3, 5683,
                                                                       5686, 22606, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47776, 3, 5686,
                                                                       5689, 22612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47786, 3, 5689,
                                                                       5692, 22618, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47796, 3, 5698,
                                                                       5701, 22624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47806, 3, 5701,
                                                                       5704, 22630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47816, 3, 5704,
                                                                       5707, 22636, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47826, 3, 5707,
                                                                       5710, 22642, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47836, 3, 5710,
                                                                       5713, 22648, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47846, 3, 5713,
                                                                       5716, 22654, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47856, 3, 5716,
                                                                       5719, 22660, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47866, 3, 5719,
                                                                       5722, 22666, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47876, 3, 5722,
                                                                       5725, 22672, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47886, 3, 5725,
                                                                       5728, 22678, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47896, 3, 5728,
                                                                       5731, 22684, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47906, 3, 5731,
                                                                       5734, 22690, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47916, 3, 5734,
                                                                       5737, 22696, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 47926, 3, 5737,
                                                                       5740, 22702, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47936, 0, 3,
                                                                       47656, 22540, 47666,
                                                                       22708, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47966, 0, 3,
                                                                       47666, 22546, 47676,
                                                                       22726, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 47996, 0, 3,
                                                                       47676, 22552, 47686,
                                                                       22744, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48026, 0, 3,
                                                                       47686, 22558, 47696,
                                                                       22762, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48056, 0, 3,
                                                                       47696, 22564, 47706,
                                                                       22780, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48086, 0, 3,
                                                                       47706, 22570, 47716,
                                                                       22798, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48116, 0, 3,
                                                                       47716, 22576, 47726,
                                                                       22816, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48146, 0, 3,
                                                                       47726, 22582, 47736,
                                                                       22834, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48176, 0, 3,
                                                                       47736, 22588, 47746,
                                                                       22852, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48206, 0, 3,
                                                                       47746, 22594, 47756,
                                                                       22870, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48236, 0, 3,
                                                                       47756, 22600, 47766,
                                                                       22888, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48266, 0, 3,
                                                                       47766, 22606, 47776,
                                                                       22906, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48296, 0, 3,
                                                                       47776, 22612, 47786,
                                                                       22924, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48326, 0, 3,
                                                                       47796, 22624, 47806,
                                                                       22942, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48356, 0, 3,
                                                                       47806, 22630, 47816,
                                                                       22960, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48386, 0, 3,
                                                                       47816, 22636, 47826,
                                                                       22978, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48416, 0, 3,
                                                                       47826, 22642, 47836,
                                                                       22996, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48446, 0, 3,
                                                                       47836, 22648, 47846,
                                                                       23014, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48476, 0, 3,
                                                                       47846, 22654, 47856,
                                                                       23032, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48506, 0, 3,
                                                                       47856, 22660, 47866,
                                                                       23050, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48536, 0, 3,
                                                                       47866, 22666, 47876,
                                                                       23068, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48566, 0, 3,
                                                                       47876, 22672, 47886,
                                                                       23086, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48596, 0, 3,
                                                                       47886, 22678, 47896,
                                                                       23104, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48626, 0, 3,
                                                                       47896, 22684, 47906,
                                                                       23122, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48656, 0, 3,
                                                                       47906, 22690, 47916,
                                                                       23140, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 48686, 0, 3,
                                                                       47916, 22696, 47926,
                                                                       23158, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48716, 0, 3,
                                                                       47936, 22708, 47966, 5980,
                                                                       5998, 23176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48776, 0, 3,
                                                                       47966, 22726, 47996, 5998,
                                                                       6016, 23212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48836, 0, 3,
                                                                       47996, 22744, 48026, 6016,
                                                                       6034, 23248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48896, 0, 3,
                                                                       48026, 22762, 48056, 6034,
                                                                       6052, 23284, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 48956, 0, 3,
                                                                       48056, 22780, 48086, 6052,
                                                                       6070, 23320, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49016, 0, 3,
                                                                       48086, 22798, 48116, 6070,
                                                                       6088, 23356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49076, 0, 3,
                                                                       48116, 22816, 48146, 6088,
                                                                       6106, 23392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49136, 0, 3,
                                                                       48146, 22834, 48176, 6106,
                                                                       6124, 23428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49196, 0, 3,
                                                                       48176, 22852, 48206, 6124,
                                                                       6142, 23464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49256, 0, 3,
                                                                       48206, 22870, 48236, 6142,
                                                                       6160, 23500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49316, 0, 3,
                                                                       48236, 22888, 48266, 6160,
                                                                       6178, 23536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49376, 0, 3,
                                                                       48266, 22906, 48296, 6178,
                                                                       6196, 23572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49436, 0, 3,
                                                                       48326, 22942, 48356, 6232,
                                                                       6250, 23608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49496, 0, 3,
                                                                       48356, 22960, 48386, 6250,
                                                                       6268, 23644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49556, 0, 3,
                                                                       48386, 22978, 48416, 6268,
                                                                       6286, 23680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49616, 0, 3,
                                                                       48416, 22996, 48446, 6286,
                                                                       6304, 23716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49676, 0, 3,
                                                                       48446, 23014, 48476, 6304,
                                                                       6322, 23752, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49736, 0, 3,
                                                                       48476, 23032, 48506, 6322,
                                                                       6340, 23788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49796, 0, 3,
                                                                       48506, 23050, 48536, 6340,
                                                                       6358, 23824, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49856, 0, 3,
                                                                       48536, 23068, 48566, 6358,
                                                                       6376, 23860, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49916, 0, 3,
                                                                       48566, 23086, 48596, 6376,
                                                                       6394, 23896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 49976, 0, 3,
                                                                       48596, 23104, 48626, 6394,
                                                                       6412, 23932, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 50036, 0, 3,
                                                                       48626, 23122, 48656, 6412,
                                                                       6430, 23968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 50096, 0, 3,
                                                                       48656, 23140, 48686, 6430,
                                                                       6448, 24004, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50156, 0, 3,
                                                                       48716, 23176, 48776, 6484,
                                                                       6514, 24040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50256, 0, 3,
                                                                       48776, 23212, 48836, 6514,
                                                                       6544, 24100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50356, 0, 3,
                                                                       48836, 23248, 48896, 6544,
                                                                       6574, 24160, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50456, 0, 3,
                                                                       48896, 23284, 48956, 6574,
                                                                       6604, 24220, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50556, 0, 3,
                                                                       48956, 23320, 49016, 6604,
                                                                       6634, 24280, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50656, 0, 3,
                                                                       49016, 23356, 49076, 6634,
                                                                       6664, 24340, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50756, 0, 3,
                                                                       49076, 23392, 49136, 6664,
                                                                       6694, 24400, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50856, 0, 3,
                                                                       49136, 23428, 49196, 6694,
                                                                       6724, 24460, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 50956, 0, 3,
                                                                       49196, 23464, 49256, 6724,
                                                                       6754, 24520, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51056, 0, 3,
                                                                       49256, 23500, 49316, 6754,
                                                                       6784, 24580, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51156, 0, 3,
                                                                       49316, 23536, 49376, 6784,
                                                                       6814, 24640, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51256, 0, 3,
                                                                       49436, 23608, 49496, 6874,
                                                                       6904, 24700, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51356, 0, 3,
                                                                       49496, 23644, 49556, 6904,
                                                                       6934, 24760, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51456, 0, 3,
                                                                       49556, 23680, 49616, 6934,
                                                                       6964, 24820, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51556, 0, 3,
                                                                       49616, 23716, 49676, 6964,
                                                                       6994, 24880, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51656, 0, 3,
                                                                       49676, 23752, 49736, 6994,
                                                                       7024, 24940, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51756, 0, 3,
                                                                       49736, 23788, 49796, 7024,
                                                                       7054, 25000, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51856, 0, 3,
                                                                       49796, 23824, 49856, 7054,
                                                                       7084, 25060, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 51956, 0, 3,
                                                                       49856, 23860, 49916, 7084,
                                                                       7114, 25120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 52056, 0, 3,
                                                                       49916, 23896, 49976, 7114,
                                                                       7144, 25180, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 52156, 0, 3,
                                                                       49976, 23932, 50036, 7144,
                                                                       7174, 25240, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 52256, 0, 3,
                                                                       50036, 23968, 50096, 7174,
                                                                       7204, 25300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 52356, 0, 3,
                                                                       50156, 24040, 50256, 7264,
                                                                       7309, 25360, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 52506, 0, 3,
                                                                       50256, 24100, 50356, 7309,
                                                                       7354, 25450, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 52656, 0, 3,
                                                                       50356, 24160, 50456, 7354,
                                                                       7399, 25540, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 52806, 0, 3,
                                                                       50456, 24220, 50556, 7399,
                                                                       7444, 25630, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 52956, 0, 3,
                                                                       50556, 24280, 50656, 7444,
                                                                       7489, 25720, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53106, 0, 3,
                                                                       50656, 24340, 50756, 7489,
                                                                       7534, 25810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53256, 0, 3,
                                                                       50756, 24400, 50856, 7534,
                                                                       7579, 25900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53406, 0, 3,
                                                                       50856, 24460, 50956, 7579,
                                                                       7624, 25990, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53556, 0, 3,
                                                                       50956, 24520, 51056, 7624,
                                                                       7669, 26080, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53706, 0, 3,
                                                                       51056, 24580, 51156, 7669,
                                                                       7714, 26170, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 53856, 0, 3,
                                                                       51256, 24700, 51356, 7804,
                                                                       7849, 26260, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54006, 0, 3,
                                                                       51356, 24760, 51456, 7849,
                                                                       7894, 26350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54156, 0, 3,
                                                                       51456, 24820, 51556, 7894,
                                                                       7939, 26440, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54306, 0, 3,
                                                                       51556, 24880, 51656, 7939,
                                                                       7984, 26530, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54456, 0, 3,
                                                                       51656, 24940, 51756, 7984,
                                                                       8029, 26620, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54606, 0, 3,
                                                                       51756, 25000, 51856, 8029,
                                                                       8074, 26710, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54756, 0, 3,
                                                                       51856, 25060, 51956, 8074,
                                                                       8119, 26800, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 54906, 0, 3,
                                                                       51956, 25120, 52056, 8119,
                                                                       8164, 26890, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 55056, 0, 3,
                                                                       52056, 25180, 52156, 8164,
                                                                       8209, 26980, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 55206, 0, 3,
                                                                       52156, 25240, 52256, 8209,
                                                                       8254, 27070, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 55356, 0, 3,
                                                                       52356, 25360, 52506, 8344,
                                                                       8407, 27160, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 55566, 0, 3,
                                                                       52506, 25450, 52656, 8407,
                                                                       8470, 27286, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 55776, 0, 3,
                                                                       52656, 25540, 52806, 8470,
                                                                       8533, 27412, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 55986, 0, 3,
                                                                       52806, 25630, 52956, 8533,
                                                                       8596, 27538, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 56196, 0, 3,
                                                                       52956, 25720, 53106, 8596,
                                                                       8659, 27664, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 56406, 0, 3,
                                                                       53106, 25810, 53256, 8659,
                                                                       8722, 27790, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 56616, 0, 3,
                                                                       53256, 25900, 53406, 8722,
                                                                       8785, 27916, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 56826, 0, 3,
                                                                       53406, 25990, 53556, 8785,
                                                                       8848, 28042, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 57036, 0, 3,
                                                                       53556, 26080, 53706, 8848,
                                                                       8911, 28168, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 57246, 0, 3,
                                                                       53856, 26260, 54006, 9037,
                                                                       9100, 28294, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 57456, 0, 3,
                                                                       54006, 26350, 54156, 9100,
                                                                       9163, 28420, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 57666, 0, 3,
                                                                       54156, 26440, 54306, 9163,
                                                                       9226, 28546, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 57876, 0, 3,
                                                                       54306, 26530, 54456, 9226,
                                                                       9289, 28672, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 58086, 0, 3,
                                                                       54456, 26620, 54606, 9289,
                                                                       9352, 28798, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 58296, 0, 3,
                                                                       54606, 26710, 54756, 9352,
                                                                       9415, 28924, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 58506, 0, 3,
                                                                       54756, 26800, 54906, 9415,
                                                                       9478, 29050, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 58716, 0, 3,
                                                                       54906, 26890, 55056, 9478,
                                                                       9541, 29176, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 58926, 0, 3,
                                                                       55056, 26980, 55206, 9541,
                                                                       9604, 29302, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 59136, 0, 3,
                                                                       55356, 27160, 55566, 9730,
                                                                       9814, 29428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 59416, 0, 3,
                                                                       55566, 27286, 55776, 9814,
                                                                       9898, 29596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 59696, 0, 3,
                                                                       55776, 27412, 55986, 9898,
                                                                       9982, 29764, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 59976, 0, 3,
                                                                       55986, 27538, 56196, 9982,
                                                                       10066, 29932, ncols,
                                                                       gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 60256, 0, 3,
                                                                       56196, 27664, 56406,
                                                                       10066, 10150, 30100,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 60536, 0, 3,
                                                                       56406, 27790, 56616,
                                                                       10150, 10234, 30268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 60816, 0, 3,
                                                                       56616, 27916, 56826,
                                                                       10234, 10318, 30436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 61096, 0, 3,
                                                                       56826, 28042, 57036,
                                                                       10318, 10402, 30604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 61376, 0, 3,
                                                                       57246, 28294, 57456,
                                                                       10570, 10654, 30772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 61656, 0, 3,
                                                                       57456, 28420, 57666,
                                                                       10654, 10738, 30940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 61936, 0, 3,
                                                                       57666, 28546, 57876,
                                                                       10738, 10822, 31108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 62216, 0, 3,
                                                                       57876, 28672, 58086,
                                                                       10822, 10906, 31276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 62496, 0, 3,
                                                                       58086, 28798, 58296,
                                                                       10906, 10990, 31444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 62776, 0, 3,
                                                                       58296, 28924, 58506,
                                                                       10990, 11074, 31612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 63056, 0, 3,
                                                                       58506, 29050, 58716,
                                                                       11074, 11158, 31780,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 63336, 0, 3,
                                                                       58716, 29176, 58926,
                                                                       11158, 11242, 31948,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 63616, 0, 3,
                                                                       59136, 29428, 59416,
                                                                       11410, 11518, 32116,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 63976, 0, 3,
                                                                       59416, 29596, 59696,
                                                                       11518, 11626, 32332,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 64336, 0, 3,
                                                                       59696, 29764, 59976,
                                                                       11626, 11734, 32548,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 64696, 0, 3,
                                                                       59976, 29932, 60256,
                                                                       11734, 11842, 32764,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 65056, 0, 3,
                                                                       60256, 30100, 60536,
                                                                       11842, 11950, 32980,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 65416, 0, 3,
                                                                       60536, 30268, 60816,
                                                                       11950, 12058, 33196,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 65776, 0, 3,
                                                                       60816, 30436, 61096,
                                                                       12058, 12166, 33412,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 66136, 0, 3,
                                                                       61376, 30772, 61656,
                                                                       12382, 12490, 33628,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 66496, 0, 3,
                                                                       61656, 30940, 61936,
                                                                       12490, 12598, 33844,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 66856, 0, 3,
                                                                       61936, 31108, 62216,
                                                                       12598, 12706, 34060,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 67216, 0, 3,
                                                                       62216, 31276, 62496,
                                                                       12706, 12814, 34276,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 67576, 0, 3,
                                                                       62496, 31444, 62776,
                                                                       12814, 12922, 34492,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 67936, 0, 3,
                                                                       62776, 31612, 63056,
                                                                       12922, 13030, 34708,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 68296, 0, 3,
                                                                       63056, 31780, 63336,
                                                                       13030, 13138, 34924,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 68656, 0, 3,
                                                                       63616, 32116, 63976,
                                                                       13354, 13489, 35140,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 69106, 0, 3,
                                                                       63976, 32332, 64336,
                                                                       13489, 13624, 35410,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 69556, 0, 3,
                                                                       64336, 32548, 64696,
                                                                       13624, 13759, 35680,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 70006, 0, 3,
                                                                       64696, 32764, 65056,
                                                                       13759, 13894, 35950,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 70456, 0, 3,
                                                                       65056, 32980, 65416,
                                                                       13894, 14029, 36220,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 70906, 0, 3,
                                                                       65416, 33196, 65776,
                                                                       14029, 14164, 36490,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 71356, 0, 3,
                                                                       66136, 33628, 66496,
                                                                       14434, 14569, 36760,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 71806, 0, 3,
                                                                       66496, 33844, 66856,
                                                                       14569, 14704, 37030,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 72256, 0, 3,
                                                                       66856, 34060, 67216,
                                                                       14704, 14839, 37300,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 72706, 0, 3,
                                                                       67216, 34276, 67576,
                                                                       14839, 14974, 37570,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 73156, 0, 3,
                                                                       67576, 34492, 67936,
                                                                       14974, 15109, 37840,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 73606, 0, 3,
                                                                       67936, 34708, 68296,
                                                                       15109, 15244, 38110,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 74056, 0, 3,
                                                                       68656, 35140, 69106,
                                                                       15514, 15679, 38380,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 74606, 0, 3,
                                                                       69106, 35410, 69556,
                                                                       15679, 15844, 38710,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 75156, 0, 3,
                                                                       69556, 35680, 70006,
                                                                       15844, 16009, 39040,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 75706, 0, 3,
                                                                       70006, 35950, 70456,
                                                                       16009, 16174, 39370,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 76256, 0, 3,
                                                                       70456, 36220, 70906,
                                                                       16174, 16339, 39700,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 76806, 0, 3,
                                                                       71356, 36760, 71806,
                                                                       16669, 16834, 40030,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 77356, 0, 3,
                                                                       71806, 37030, 72256,
                                                                       16834, 16999, 40360,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 77906, 0, 3,
                                                                       72256, 37300, 72706,
                                                                       16999, 17164, 40690,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 78456, 0, 3,
                                                                       72706, 37570, 73156,
                                                                       17164, 17329, 41020,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 79006, 0, 3,
                                                                       73156, 37840, 73606,
                                                                       17329, 17494, 41350,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 79556, 0, 3,
                                                                       74056, 38380, 74606,
                                                                       17824, 18022, 41680,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 80216, 0, 3,
                                                                       74606, 38710, 75156,
                                                                       18022, 18220, 42076,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 80876, 0, 3,
                                                                       75156, 39040, 75706,
                                                                       18220, 18418, 42472,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 81536, 0, 3,
                                                                       75706, 39370, 76256,
                                                                       18418, 18616, 42868,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 82196, 0, 3,
                                                                       76806, 40030, 77356,
                                                                       19012, 19210, 43264,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 82856, 0, 3,
                                                                       77356, 40360, 77906,
                                                                       19210, 19408, 43660,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 83516, 0, 3,
                                                                       77906, 40690, 78456,
                                                                       19408, 19606, 44056,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 84176, 0, 3,
                                                                       78456, 41020, 79006,
                                                                       19606, 19804, 44452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 84836, 0, 3,
                                                                       79556, 41680, 80216,
                                                                       20200, 20434, 44848,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 85616, 0, 3,
                                                                       80216, 42076, 80876,
                                                                       20434, 20668, 45316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 86396, 0, 3,
                                                                       80876, 42472, 81536,
                                                                       20668, 20902, 45784,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 87176, 0, 3,
                                                                       82196, 43264, 82856,
                                                                       21370, 21604, 46252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 87956, 0, 3,
                                                                       82856, 43660, 83516,
                                                                       21604, 21838, 46720,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 88736, 0, 3,
                                                                       83516, 44056, 84176,
                                                                       21838, 22072, 47188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89516, 3, 22540,
                                                                       22546, 47676, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89531, 3, 22546,
                                                                       22552, 47686, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89546, 3, 22552,
                                                                       22558, 47696, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89561, 3, 22558,
                                                                       22564, 47706, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89576, 3, 22564,
                                                                       22570, 47716, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89591, 3, 22570,
                                                                       22576, 47726, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89606, 3, 22576,
                                                                       22582, 47736, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89621, 3, 22582,
                                                                       22588, 47746, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89636, 3, 22588,
                                                                       22594, 47756, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89651, 3, 22594,
                                                                       22600, 47766, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89666, 3, 22600,
                                                                       22606, 47776, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89681, 3, 22606,
                                                                       22612, 47786, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89696, 3, 22624,
                                                                       22630, 47816, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89711, 3, 22630,
                                                                       22636, 47826, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89726, 3, 22636,
                                                                       22642, 47836, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89741, 3, 22642,
                                                                       22648, 47846, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89756, 3, 22648,
                                                                       22654, 47856, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89771, 3, 22654,
                                                                       22660, 47866, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89786, 3, 22660,
                                                                       22666, 47876, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89801, 3, 22666,
                                                                       22672, 47886, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89816, 3, 22672,
                                                                       22678, 47896, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89831, 3, 22678,
                                                                       22684, 47906, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89846, 3, 22684,
                                                                       22690, 47916, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 89861, 3, 22690,
                                                                       22696, 47926, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 89876, 0, 3,
                                                                       89516, 47676, 89531,
                                                                       22708, 22726, 47996,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 89921, 0, 3,
                                                                       89531, 47686, 89546,
                                                                       22726, 22744, 48026,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 89966, 0, 3,
                                                                       89546, 47696, 89561,
                                                                       22744, 22762, 48056,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90011, 0, 3,
                                                                       89561, 47706, 89576,
                                                                       22762, 22780, 48086,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90056, 0, 3,
                                                                       89576, 47716, 89591,
                                                                       22780, 22798, 48116,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90101, 0, 3,
                                                                       89591, 47726, 89606,
                                                                       22798, 22816, 48146,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90146, 0, 3,
                                                                       89606, 47736, 89621,
                                                                       22816, 22834, 48176,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90191, 0, 3,
                                                                       89621, 47746, 89636,
                                                                       22834, 22852, 48206,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90236, 0, 3,
                                                                       89636, 47756, 89651,
                                                                       22852, 22870, 48236,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90281, 0, 3,
                                                                       89651, 47766, 89666,
                                                                       22870, 22888, 48266,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90326, 0, 3,
                                                                       89666, 47776, 89681,
                                                                       22888, 22906, 48296,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90371, 0, 3,
                                                                       89696, 47816, 89711,
                                                                       22942, 22960, 48386,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90416, 0, 3,
                                                                       89711, 47826, 89726,
                                                                       22960, 22978, 48416,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90461, 0, 3,
                                                                       89726, 47836, 89741,
                                                                       22978, 22996, 48446,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90506, 0, 3,
                                                                       89741, 47846, 89756,
                                                                       22996, 23014, 48476,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90551, 0, 3,
                                                                       89756, 47856, 89771,
                                                                       23014, 23032, 48506,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90596, 0, 3,
                                                                       89771, 47866, 89786,
                                                                       23032, 23050, 48536,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90641, 0, 3,
                                                                       89786, 47876, 89801,
                                                                       23050, 23068, 48566,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90686, 0, 3,
                                                                       89801, 47886, 89816,
                                                                       23068, 23086, 48596,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90731, 0, 3,
                                                                       89816, 47896, 89831,
                                                                       23086, 23104, 48626,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90776, 0, 3,
                                                                       89831, 47906, 89846,
                                                                       23104, 23122, 48656,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 90821, 0, 3,
                                                                       89846, 47916, 89861,
                                                                       23122, 23140, 48686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 90866, 0, 3,
                                                                       89876, 47996, 89921,
                                                                       23176, 23212, 48836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 90956, 0, 3,
                                                                       89921, 48026, 89966,
                                                                       23212, 23248, 48896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91046, 0, 3,
                                                                       89966, 48056, 90011,
                                                                       23248, 23284, 48956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91136, 0, 3,
                                                                       90011, 48086, 90056,
                                                                       23284, 23320, 49016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91226, 0, 3,
                                                                       90056, 48116, 90101,
                                                                       23320, 23356, 49076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91316, 0, 3,
                                                                       90101, 48146, 90146,
                                                                       23356, 23392, 49136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91406, 0, 3,
                                                                       90146, 48176, 90191,
                                                                       23392, 23428, 49196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91496, 0, 3,
                                                                       90191, 48206, 90236,
                                                                       23428, 23464, 49256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91586, 0, 3,
                                                                       90236, 48236, 90281,
                                                                       23464, 23500, 49316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91676, 0, 3,
                                                                       90281, 48266, 90326,
                                                                       23500, 23536, 49376,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91766, 0, 3,
                                                                       90371, 48386, 90416,
                                                                       23608, 23644, 49556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91856, 0, 3,
                                                                       90416, 48416, 90461,
                                                                       23644, 23680, 49616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 91946, 0, 3,
                                                                       90461, 48446, 90506,
                                                                       23680, 23716, 49676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92036, 0, 3,
                                                                       90506, 48476, 90551,
                                                                       23716, 23752, 49736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92126, 0, 3,
                                                                       90551, 48506, 90596,
                                                                       23752, 23788, 49796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92216, 0, 3,
                                                                       90596, 48536, 90641,
                                                                       23788, 23824, 49856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92306, 0, 3,
                                                                       90641, 48566, 90686,
                                                                       23824, 23860, 49916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92396, 0, 3,
                                                                       90686, 48596, 90731,
                                                                       23860, 23896, 49976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92486, 0, 3,
                                                                       90731, 48626, 90776,
                                                                       23896, 23932, 50036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 92576, 0, 3,
                                                                       90776, 48656, 90821,
                                                                       23932, 23968, 50096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 92666, 0, 3,
                                                                       90866, 48836, 90956,
                                                                       24040, 24100, 50356,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 92816, 0, 3,
                                                                       90956, 48896, 91046,
                                                                       24100, 24160, 50456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 92966, 0, 3,
                                                                       91046, 48956, 91136,
                                                                       24160, 24220, 50556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 93116, 0, 3,
                                                                       91136, 49016, 91226,
                                                                       24220, 24280, 50656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 93266, 0, 3,
                                                                       91226, 49076, 91316,
                                                                       24280, 24340, 50756,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 93416, 0, 3,
                                                                       91316, 49136, 91406,
                                                                       24340, 24400, 50856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 93566, 0, 3,
                                                                       91406, 49196, 91496,
                                                                       24400, 24460, 50956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 93716, 0, 3,
                                                                       91496, 49256, 91586,
                                                                       24460, 24520, 51056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 93866, 0, 3,
                                                                       91586, 49316, 91676,
                                                                       24520, 24580, 51156,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94016, 0, 3,
                                                                       91766, 49556, 91856,
                                                                       24700, 24760, 51456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94166, 0, 3,
                                                                       91856, 49616, 91946,
                                                                       24760, 24820, 51556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94316, 0, 3,
                                                                       91946, 49676, 92036,
                                                                       24820, 24880, 51656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94466, 0, 3,
                                                                       92036, 49736, 92126,
                                                                       24880, 24940, 51756,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94616, 0, 3,
                                                                       92126, 49796, 92216,
                                                                       24940, 25000, 51856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94766, 0, 3,
                                                                       92216, 49856, 92306,
                                                                       25000, 25060, 51956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 94916, 0, 3,
                                                                       92306, 49916, 92396,
                                                                       25060, 25120, 52056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 95066, 0, 3,
                                                                       92396, 49976, 92486,
                                                                       25120, 25180, 52156,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 95216, 0, 3,
                                                                       92486, 50036, 92576,
                                                                       25180, 25240, 52256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 95366, 0, 3,
                                                                       92666, 50356, 92816,
                                                                       25360, 25450, 52656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 95591, 0, 3,
                                                                       92816, 50456, 92966,
                                                                       25450, 25540, 52806,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 95816, 0, 3,
                                                                       92966, 50556, 93116,
                                                                       25540, 25630, 52956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 96041, 0, 3,
                                                                       93116, 50656, 93266,
                                                                       25630, 25720, 53106,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 96266, 0, 3,
                                                                       93266, 50756, 93416,
                                                                       25720, 25810, 53256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 96491, 0, 3,
                                                                       93416, 50856, 93566,
                                                                       25810, 25900, 53406,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 96716, 0, 3,
                                                                       93566, 50956, 93716,
                                                                       25900, 25990, 53556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 96941, 0, 3,
                                                                       93716, 51056, 93866,
                                                                       25990, 26080, 53706,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 97166, 0, 3,
                                                                       94016, 51456, 94166,
                                                                       26260, 26350, 54156,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 97391, 0, 3,
                                                                       94166, 51556, 94316,
                                                                       26350, 26440, 54306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 97616, 0, 3,
                                                                       94316, 51656, 94466,
                                                                       26440, 26530, 54456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 97841, 0, 3,
                                                                       94466, 51756, 94616,
                                                                       26530, 26620, 54606,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 98066, 0, 3,
                                                                       94616, 51856, 94766,
                                                                       26620, 26710, 54756,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 98291, 0, 3,
                                                                       94766, 51956, 94916,
                                                                       26710, 26800, 54906,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 98516, 0, 3,
                                                                       94916, 52056, 95066,
                                                                       26800, 26890, 55056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 98741, 0, 3,
                                                                       95066, 52156, 95216,
                                                                       26890, 26980, 55206,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 98966, 0, 3,
                                                                       95366, 52656, 95591,
                                                                       27160, 27286, 55776,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 99281, 0, 3,
                                                                       95591, 52806, 95816,
                                                                       27286, 27412, 55986,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 99596, 0, 3,
                                                                       95816, 52956, 96041,
                                                                       27412, 27538, 56196,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 99911, 0, 3,
                                                                       96041, 53106, 96266,
                                                                       27538, 27664, 56406,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 100226, 0, 3,
                                                                       96266, 53256, 96491,
                                                                       27664, 27790, 56616,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 100541, 0, 3,
                                                                       96491, 53406, 96716,
                                                                       27790, 27916, 56826,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 100856, 0, 3,
                                                                       96716, 53556, 96941,
                                                                       27916, 28042, 57036,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 101171, 0, 3,
                                                                       97166, 54156, 97391,
                                                                       28294, 28420, 57666,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 101486, 0, 3,
                                                                       97391, 54306, 97616,
                                                                       28420, 28546, 57876,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 101801, 0, 3,
                                                                       97616, 54456, 97841,
                                                                       28546, 28672, 58086,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 102116, 0, 3,
                                                                       97841, 54606, 98066,
                                                                       28672, 28798, 58296,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 102431, 0, 3,
                                                                       98066, 54756, 98291,
                                                                       28798, 28924, 58506,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 102746, 0, 3,
                                                                       98291, 54906, 98516,
                                                                       28924, 29050, 58716,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 103061, 0, 3,
                                                                       98516, 55056, 98741,
                                                                       29050, 29176, 58926,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 103376, 0, 3,
                                                                       98966, 55776, 99281,
                                                                       29428, 29596, 59696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 103796, 0, 3,
                                                                       99281, 55986, 99596,
                                                                       29596, 29764, 59976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 104216, 0, 3,
                                                                       99596, 56196, 99911,
                                                                       29764, 29932, 60256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 104636, 0, 3,
                                                                       99911, 56406, 100226,
                                                                       29932, 30100, 60536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 105056, 0, 3,
                                                                       100226, 56616, 100541,
                                                                       30100, 30268, 60816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 105476, 0, 3,
                                                                       100541, 56826, 100856,
                                                                       30268, 30436, 61096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 105896, 0, 3,
                                                                       101171, 57666, 101486,
                                                                       30772, 30940, 61936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 106316, 0, 3,
                                                                       101486, 57876, 101801,
                                                                       30940, 31108, 62216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 106736, 0, 3,
                                                                       101801, 58086, 102116,
                                                                       31108, 31276, 62496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 107156, 0, 3,
                                                                       102116, 58296, 102431,
                                                                       31276, 31444, 62776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 107576, 0, 3,
                                                                       102431, 58506, 102746,
                                                                       31444, 31612, 63056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 107996, 0, 3,
                                                                       102746, 58716, 103061,
                                                                       31612, 31780, 63336,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 108416, 0, 3,
                                                                       103376, 59696, 103796,
                                                                       32116, 32332, 64336,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 108956, 0, 3,
                                                                       103796, 59976, 104216,
                                                                       32332, 32548, 64696,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 109496, 0, 3,
                                                                       104216, 60256, 104636,
                                                                       32548, 32764, 65056,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 110036, 0, 3,
                                                                       104636, 60536, 105056,
                                                                       32764, 32980, 65416,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 110576, 0, 3,
                                                                       105056, 60816, 105476,
                                                                       32980, 33196, 65776,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 111116, 0, 3,
                                                                       105896, 61936, 106316,
                                                                       33628, 33844, 66856,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 111656, 0, 3,
                                                                       106316, 62216, 106736,
                                                                       33844, 34060, 67216,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 112196, 0, 3,
                                                                       106736, 62496, 107156,
                                                                       34060, 34276, 67576,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 112736, 0, 3,
                                                                       107156, 62776, 107576,
                                                                       34276, 34492, 67936,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 113276, 0, 3,
                                                                       107576, 63056, 107996,
                                                                       34492, 34708, 68296,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 113816, 0, 3,
                                                                       108416, 64336, 108956,
                                                                       35140, 35410, 69556,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 114491, 0, 3,
                                                                       108956, 64696, 109496,
                                                                       35410, 35680, 70006,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 115166, 0, 3,
                                                                       109496, 65056, 110036,
                                                                       35680, 35950, 70456,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 115841, 0, 3,
                                                                       110036, 65416, 110576,
                                                                       35950, 36220, 70906,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 116516, 0, 3,
                                                                       111116, 66856, 111656,
                                                                       36760, 37030, 72256,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 117191, 0, 3,
                                                                       111656, 67216, 112196,
                                                                       37030, 37300, 72706,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 117866, 0, 3,
                                                                       112196, 67576, 112736,
                                                                       37300, 37570, 73156,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 118541, 0, 3,
                                                                       112736, 67936, 113276,
                                                                       37570, 37840, 73606,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 119216, 0, 3,
                                                                       113816, 69556, 114491,
                                                                       38380, 38710, 75156,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 120041, 0, 3,
                                                                       114491, 70006, 115166,
                                                                       38710, 39040, 75706,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 120866, 0, 3,
                                                                       115166, 70456, 115841,
                                                                       39040, 39370, 76256,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 121691, 0, 3,
                                                                       116516, 72256, 117191,
                                                                       40030, 40360, 77906,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 122516, 0, 3,
                                                                       117191, 72706, 117866,
                                                                       40360, 40690, 78456,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 123341, 0, 3,
                                                                       117866, 73156, 118541,
                                                                       40690, 41020, 79006,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 124166, 0, 3,
                                                                       119216, 75156, 120041,
                                                                       41680, 42076, 80876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 125156, 0, 3,
                                                                       120041, 75706, 120866,
                                                                       42076, 42472, 81536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 126146, 0, 3,
                                                                       121691, 77906, 122516,
                                                                       43264, 43660, 83516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 127136, 0, 3,
                                                                       122516, 78456, 123341,
                                                                       43660, 44056, 84176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 128126, 0, 3,
                                                                       124166, 80876, 125156,
                                                                       44848, 45316, 86396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 129296, 0, 3,
                                                                       126146, 83516, 127136,
                                                                       46252, 46720, 88736,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130466, 3, 47656,
                                                                       47666, 89516, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130487, 3, 47666,
                                                                       47676, 89531, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130508, 3, 47676,
                                                                       47686, 89546, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130529, 3, 47686,
                                                                       47696, 89561, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130550, 3, 47696,
                                                                       47706, 89576, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130571, 3, 47706,
                                                                       47716, 89591, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130592, 3, 47716,
                                                                       47726, 89606, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130613, 3, 47726,
                                                                       47736, 89621, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130634, 3, 47736,
                                                                       47746, 89636, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130655, 3, 47746,
                                                                       47756, 89651, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130676, 3, 47756,
                                                                       47766, 89666, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130697, 3, 47766,
                                                                       47776, 89681, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130718, 3, 47796,
                                                                       47806, 89696, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130739, 3, 47806,
                                                                       47816, 89711, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130760, 3, 47816,
                                                                       47826, 89726, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130781, 3, 47826,
                                                                       47836, 89741, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130802, 3, 47836,
                                                                       47846, 89756, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130823, 3, 47846,
                                                                       47856, 89771, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130844, 3, 47856,
                                                                       47866, 89786, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130865, 3, 47866,
                                                                       47876, 89801, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130886, 3, 47876,
                                                                       47886, 89816, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130907, 3, 47886,
                                                                       47896, 89831, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130928, 3, 47896,
                                                                       47906, 89846, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 130949, 3, 47906,
                                                                       47916, 89861, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 130970, 0, 3,
                                                                       130466, 89516, 130487,
                                                                       47936, 47966, 89876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131033, 0, 3,
                                                                       130487, 89531, 130508,
                                                                       47966, 47996, 89921,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131096, 0, 3,
                                                                       130508, 89546, 130529,
                                                                       47996, 48026, 89966,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131159, 0, 3,
                                                                       130529, 89561, 130550,
                                                                       48026, 48056, 90011,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131222, 0, 3,
                                                                       130550, 89576, 130571,
                                                                       48056, 48086, 90056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131285, 0, 3,
                                                                       130571, 89591, 130592,
                                                                       48086, 48116, 90101,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131348, 0, 3,
                                                                       130592, 89606, 130613,
                                                                       48116, 48146, 90146,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131411, 0, 3,
                                                                       130613, 89621, 130634,
                                                                       48146, 48176, 90191,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131474, 0, 3,
                                                                       130634, 89636, 130655,
                                                                       48176, 48206, 90236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131537, 0, 3,
                                                                       130655, 89651, 130676,
                                                                       48206, 48236, 90281,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131600, 0, 3,
                                                                       130676, 89666, 130697,
                                                                       48236, 48266, 90326,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131663, 0, 3,
                                                                       130718, 89696, 130739,
                                                                       48326, 48356, 90371,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131726, 0, 3,
                                                                       130739, 89711, 130760,
                                                                       48356, 48386, 90416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131789, 0, 3,
                                                                       130760, 89726, 130781,
                                                                       48386, 48416, 90461,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131852, 0, 3,
                                                                       130781, 89741, 130802,
                                                                       48416, 48446, 90506,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131915, 0, 3,
                                                                       130802, 89756, 130823,
                                                                       48446, 48476, 90551,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 131978, 0, 3,
                                                                       130823, 89771, 130844,
                                                                       48476, 48506, 90596,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 132041, 0, 3,
                                                                       130844, 89786, 130865,
                                                                       48506, 48536, 90641,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 132104, 0, 3,
                                                                       130865, 89801, 130886,
                                                                       48536, 48566, 90686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 132167, 0, 3,
                                                                       130886, 89816, 130907,
                                                                       48566, 48596, 90731,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 132230, 0, 3,
                                                                       130907, 89831, 130928,
                                                                       48596, 48626, 90776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 132293, 0, 3,
                                                                       130928, 89846, 130949,
                                                                       48626, 48656, 90821,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 132356, 0, 3,
                                                                       130970, 89876, 131033,
                                                                       48716, 48776, 90866,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 132482, 0, 3,
                                                                       131033, 89921, 131096,
                                                                       48776, 48836, 90956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 132608, 0, 3,
                                                                       131096, 89966, 131159,
                                                                       48836, 48896, 91046,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 132734, 0, 3,
                                                                       131159, 90011, 131222,
                                                                       48896, 48956, 91136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 132860, 0, 3,
                                                                       131222, 90056, 131285,
                                                                       48956, 49016, 91226,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 132986, 0, 3,
                                                                       131285, 90101, 131348,
                                                                       49016, 49076, 91316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 133112, 0, 3,
                                                                       131348, 90146, 131411,
                                                                       49076, 49136, 91406,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 133238, 0, 3,
                                                                       131411, 90191, 131474,
                                                                       49136, 49196, 91496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 133364, 0, 3,
                                                                       131474, 90236, 131537,
                                                                       49196, 49256, 91586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 133490, 0, 3,
                                                                       131537, 90281, 131600,
                                                                       49256, 49316, 91676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 133616, 0, 3,
                                                                       131663, 90371, 131726,
                                                                       49436, 49496, 91766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 133742, 0, 3,
                                                                       131726, 90416, 131789,
                                                                       49496, 49556, 91856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 133868, 0, 3,
                                                                       131789, 90461, 131852,
                                                                       49556, 49616, 91946,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 133994, 0, 3,
                                                                       131852, 90506, 131915,
                                                                       49616, 49676, 92036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 134120, 0, 3,
                                                                       131915, 90551, 131978,
                                                                       49676, 49736, 92126,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 134246, 0, 3,
                                                                       131978, 90596, 132041,
                                                                       49736, 49796, 92216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 134372, 0, 3,
                                                                       132041, 90641, 132104,
                                                                       49796, 49856, 92306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 134498, 0, 3,
                                                                       132104, 90686, 132167,
                                                                       49856, 49916, 92396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 134624, 0, 3,
                                                                       132167, 90731, 132230,
                                                                       49916, 49976, 92486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 134750, 0, 3,
                                                                       132230, 90776, 132293,
                                                                       49976, 50036, 92576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 134876, 0, 3,
                                                                       132356, 90866, 132482,
                                                                       50156, 50256, 92666,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 135086, 0, 3,
                                                                       132482, 90956, 132608,
                                                                       50256, 50356, 92816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 135296, 0, 3,
                                                                       132608, 91046, 132734,
                                                                       50356, 50456, 92966,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 135506, 0, 3,
                                                                       132734, 91136, 132860,
                                                                       50456, 50556, 93116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 135716, 0, 3,
                                                                       132860, 91226, 132986,
                                                                       50556, 50656, 93266,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 135926, 0, 3,
                                                                       132986, 91316, 133112,
                                                                       50656, 50756, 93416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 136136, 0, 3,
                                                                       133112, 91406, 133238,
                                                                       50756, 50856, 93566,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 136346, 0, 3,
                                                                       133238, 91496, 133364,
                                                                       50856, 50956, 93716,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 136556, 0, 3,
                                                                       133364, 91586, 133490,
                                                                       50956, 51056, 93866,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 136766, 0, 3,
                                                                       133616, 91766, 133742,
                                                                       51256, 51356, 94016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 136976, 0, 3,
                                                                       133742, 91856, 133868,
                                                                       51356, 51456, 94166,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 137186, 0, 3,
                                                                       133868, 91946, 133994,
                                                                       51456, 51556, 94316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 137396, 0, 3,
                                                                       133994, 92036, 134120,
                                                                       51556, 51656, 94466,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 137606, 0, 3,
                                                                       134120, 92126, 134246,
                                                                       51656, 51756, 94616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 137816, 0, 3,
                                                                       134246, 92216, 134372,
                                                                       51756, 51856, 94766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 138026, 0, 3,
                                                                       134372, 92306, 134498,
                                                                       51856, 51956, 94916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 138236, 0, 3,
                                                                       134498, 92396, 134624,
                                                                       51956, 52056, 95066,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 138446, 0, 3,
                                                                       134624, 92486, 134750,
                                                                       52056, 52156, 95216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 138656, 0, 3,
                                                                       134876, 92666, 135086,
                                                                       52356, 52506, 95366,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 138971, 0, 3,
                                                                       135086, 92816, 135296,
                                                                       52506, 52656, 95591,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 139286, 0, 3,
                                                                       135296, 92966, 135506,
                                                                       52656, 52806, 95816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 139601, 0, 3,
                                                                       135506, 93116, 135716,
                                                                       52806, 52956, 96041,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 139916, 0, 3,
                                                                       135716, 93266, 135926,
                                                                       52956, 53106, 96266,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 140231, 0, 3,
                                                                       135926, 93416, 136136,
                                                                       53106, 53256, 96491,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 140546, 0, 3,
                                                                       136136, 93566, 136346,
                                                                       53256, 53406, 96716,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 140861, 0, 3,
                                                                       136346, 93716, 136556,
                                                                       53406, 53556, 96941,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 141176, 0, 3,
                                                                       136766, 94016, 136976,
                                                                       53856, 54006, 97166,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 141491, 0, 3,
                                                                       136976, 94166, 137186,
                                                                       54006, 54156, 97391,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 141806, 0, 3,
                                                                       137186, 94316, 137396,
                                                                       54156, 54306, 97616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 142121, 0, 3,
                                                                       137396, 94466, 137606,
                                                                       54306, 54456, 97841,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 142436, 0, 3,
                                                                       137606, 94616, 137816,
                                                                       54456, 54606, 98066,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 142751, 0, 3,
                                                                       137816, 94766, 138026,
                                                                       54606, 54756, 98291,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 143066, 0, 3,
                                                                       138026, 94916, 138236,
                                                                       54756, 54906, 98516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 143381, 0, 3,
                                                                       138236, 95066, 138446,
                                                                       54906, 55056, 98741,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 143696, 0, 3,
                                                                       138656, 95366, 138971,
                                                                       55356, 55566, 98966,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 144137, 0, 3,
                                                                       138971, 95591, 139286,
                                                                       55566, 55776, 99281,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 144578, 0, 3,
                                                                       139286, 95816, 139601,
                                                                       55776, 55986, 99596,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 145019, 0, 3,
                                                                       139601, 96041, 139916,
                                                                       55986, 56196, 99911,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 145460, 0, 3,
                                                                       139916, 96266, 140231,
                                                                       56196, 56406, 100226,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 145901, 0, 3,
                                                                       140231, 96491, 140546,
                                                                       56406, 56616, 100541,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 146342, 0, 3,
                                                                       140546, 96716, 140861,
                                                                       56616, 56826, 100856,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 146783, 0, 3,
                                                                       141176, 97166, 141491,
                                                                       57246, 57456, 101171,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 147224, 0, 3,
                                                                       141491, 97391, 141806,
                                                                       57456, 57666, 101486,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 147665, 0, 3,
                                                                       141806, 97616, 142121,
                                                                       57666, 57876, 101801,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 148106, 0, 3,
                                                                       142121, 97841, 142436,
                                                                       57876, 58086, 102116,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 148547, 0, 3,
                                                                       142436, 98066, 142751,
                                                                       58086, 58296, 102431,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 148988, 0, 3,
                                                                       142751, 98291, 143066,
                                                                       58296, 58506, 102746,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 149429, 0, 3,
                                                                       143066, 98516, 143381,
                                                                       58506, 58716, 103061,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 149870, 0, 3,
                                                                       143696, 98966, 144137,
                                                                       59136, 59416, 103376,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 150458, 0, 3,
                                                                       144137, 99281, 144578,
                                                                       59416, 59696, 103796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 151046, 0, 3,
                                                                       144578, 99596, 145019,
                                                                       59696, 59976, 104216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 151634, 0, 3,
                                                                       145019, 99911, 145460,
                                                                       59976, 60256, 104636,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 152222, 0, 3,
                                                                       145460, 100226, 145901,
                                                                       60256, 60536, 105056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 152810, 0, 3,
                                                                       145901, 100541, 146342,
                                                                       60536, 60816, 105476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 153398, 0, 3,
                                                                       146783, 101171, 147224,
                                                                       61376, 61656, 105896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 153986, 0, 3,
                                                                       147224, 101486, 147665,
                                                                       61656, 61936, 106316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 154574, 0, 3,
                                                                       147665, 101801, 148106,
                                                                       61936, 62216, 106736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 155162, 0, 3,
                                                                       148106, 102116, 148547,
                                                                       62216, 62496, 107156,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 155750, 0, 3,
                                                                       148547, 102431, 148988,
                                                                       62496, 62776, 107576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 156338, 0, 3,
                                                                       148988, 102746, 149429,
                                                                       62776, 63056, 107996,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 156926, 0, 3,
                                                                       149870, 103376, 150458,
                                                                       63616, 63976, 108416,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 157682, 0, 3,
                                                                       150458, 103796, 151046,
                                                                       63976, 64336, 108956,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 158438, 0, 3,
                                                                       151046, 104216, 151634,
                                                                       64336, 64696, 109496,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 159194, 0, 3,
                                                                       151634, 104636, 152222,
                                                                       64696, 65056, 110036,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 159950, 0, 3,
                                                                       152222, 105056, 152810,
                                                                       65056, 65416, 110576,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 160706, 0, 3,
                                                                       153398, 105896, 153986,
                                                                       66136, 66496, 111116,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 161462, 0, 3,
                                                                       153986, 106316, 154574,
                                                                       66496, 66856, 111656,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 162218, 0, 3,
                                                                       154574, 106736, 155162,
                                                                       66856, 67216, 112196,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 162974, 0, 3,
                                                                       155162, 107156, 155750,
                                                                       67216, 67576, 112736,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 163730, 0, 3,
                                                                       155750, 107576, 156338,
                                                                       67576, 67936, 113276,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 164486, 0, 3,
                                                                       156926, 108416, 157682,
                                                                       68656, 69106, 113816,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 165431, 0, 3,
                                                                       157682, 108956, 158438,
                                                                       69106, 69556, 114491,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 166376, 0, 3,
                                                                       158438, 109496, 159194,
                                                                       69556, 70006, 115166,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 167321, 0, 3,
                                                                       159194, 110036, 159950,
                                                                       70006, 70456, 115841,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 168266, 0, 3,
                                                                       160706, 111116, 161462,
                                                                       71356, 71806, 116516,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 169211, 0, 3,
                                                                       161462, 111656, 162218,
                                                                       71806, 72256, 117191,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 170156, 0, 3,
                                                                       162218, 112196, 162974,
                                                                       72256, 72706, 117866,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 171101, 0, 3,
                                                                       162974, 112736, 163730,
                                                                       72706, 73156, 118541,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 172046, 0, 3,
                                                                       164486, 113816, 165431,
                                                                       74056, 74606, 119216,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 173201, 0, 3,
                                                                       165431, 114491, 166376,
                                                                       74606, 75156, 120041,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 174356, 0, 3,
                                                                       166376, 115166, 167321,
                                                                       75156, 75706, 120866,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 175511, 0, 3,
                                                                       168266, 116516, 169211,
                                                                       76806, 77356, 121691,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 176666, 0, 3,
                                                                       169211, 117191, 170156,
                                                                       77356, 77906, 122516,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 177821, 0, 3,
                                                                       170156, 117866, 171101,
                                                                       77906, 78456, 123341,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 178976, 0, 3,
                                                                       172046, 119216, 173201,
                                                                       79556, 80216, 124166,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 180362, 0, 3,
                                                                       173201, 120041, 174356,
                                                                       80216, 80876, 125156,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 181748, 0, 3,
                                                                       175511, 121691, 176666,
                                                                       82196, 82856, 126146,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 183134, 0, 3,
                                                                       176666, 122516, 177821,
                                                                       82856, 83516, 127136,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 184520, 0, 3,
                                                                       178976, 124166, 180362,
                                                                       84836, 85616, 128126,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 186158, 0, 3,
                                                                       181748, 126146, 183134,
                                                                       87176, 87956, 129296,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 187796, 149870, 588, ncols);

                    simdfunc::contract_primitives(buffer, 188692, 153398, 588, ncols);

                    simdfunc::contract_primitives(buffer, 189588, 156926, 756, ncols);

                    simdfunc::contract_primitives(buffer, 190740, 160706, 756, ncols);

                    simdfunc::contract_primitives(buffer, 191892, 164486, 945, ncols);

                    simdfunc::contract_primitives(buffer, 193332, 168266, 945, ncols);

                    simdfunc::contract_primitives(buffer, 194772, 172046, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 196532, 175511, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 198292, 178976, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 200404, 181748, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 202516, 184520, 1638, ncols);

                    simdfunc::contract_primitives(buffer, 205012, 186158, 1638, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 188384, 187796, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 189280, 188692, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 190344, 189588, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 191496, 190740, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 192837, 191892, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 194277, 193332, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 195927, 194772, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 197687, 196532, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 199678, 198292, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 201790, 200404, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 204154, 202516, 78, 1, nmax);

        simdtrf::transform_h_inner(buffer, 206650, 205012, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 207508, 188384, 190344, 11, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 208432, 189280, 191496, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 209356, 190344, 192837, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 210544, 191496, 194277, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 211732, 192837, 195927, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 213217, 194277, 197687, 11, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 214702, 195927, 199678, 11, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 216517, 197687, 201790, 11, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 218332, 199678, 204154, 11, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 220510, 201790, 206650, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 222688, 207508, 209356, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 224536, 208432, 210544, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 226384, 209356, 211732, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 228760, 210544, 213217, 11, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 231136, 211732, 214702, 11, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 234106, 213217, 216517, 11, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 237076, 214702, 218332, 11, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 240706, 216517, 220510, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 244336, 222688, 226384, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 247416, 224536, 228760, 11, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 250496, 226384, 231136, 11, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 254456, 228760, 234106, 11, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 258416, 231136, 237076, 11, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 263366, 234106, 240706, 11, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 268316, 244336, 250496, 11, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 272936, 247416, 254456, 11, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 277556, 250496, 258416, 11, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 283496, 254456, 263366, 11, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 289436, 268316, 277556, 11, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 295904, 272936, 283496, 11, nmax);

        simdtrf::transform_i_inner(buffer, 302372, 295904, 21, 11, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 302372, 143, nmax);

        simdtrf::transform_i_inner(buffer, 302372, 289436, 21, 11, nmax);

        simdtrf::transform_h_outer(values + 1573 * nvalues + n * npairs, nvalues, buffer, 302372,
                                   143, nmax);
    }

    for (size_t m = 0; m < 3146; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
