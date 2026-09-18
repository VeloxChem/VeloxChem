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


#include "SimdThreeCenterElectronRepulsionRsRecHIK.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"
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
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_hik_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_hik_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 608923, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 4290 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 608923, 444052, 30246, dimensions);

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
                                                            15, 16, 17, 18}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 25, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17, 18}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 77, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 80, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 83, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 86, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 89, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 92, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 95, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 98, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 101, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 104, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 107, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 110, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 113, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 116, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 119, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 122, 0, 3, 35, 36,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 125, 0, 3, 36, 37,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 128, 0, 3, 37, 38,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 131, 0, 3, 38, 39,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 134, 0, 3, 39, 40,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 137, 0, 3, 40, 41,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 140, 0, 3, 41, 42,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 143, 0, 3, 42, 43,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 146, 0, 3, 7, 8,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 152, 0, 3, 8, 9,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 158, 0, 3, 9, 10,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 164, 0, 3, 10, 11,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 170, 0, 3, 11, 12,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 176, 0, 3, 12, 13,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 182, 0, 3, 13, 14,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 188, 0, 3, 14, 15,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 194, 0, 3, 15, 16,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 200, 0, 3, 16, 17,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 206, 0, 3, 17, 18,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 212, 0, 3, 18, 19,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 218, 0, 3, 19, 20,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 224, 0, 3, 20, 21,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 230, 0, 3, 21, 22,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 236, 0, 3, 22, 23,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 242, 0, 3, 26, 27,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 248, 0, 3, 27, 28,
                                                                       98, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 254, 0, 3, 28, 29,
                                                                       101, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 260, 0, 3, 29, 30,
                                                                       104, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 266, 0, 3, 30, 31,
                                                                       107, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 272, 0, 3, 31, 32,
                                                                       110, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 278, 0, 3, 32, 33,
                                                                       113, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 284, 0, 3, 33, 34,
                                                                       116, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 290, 0, 3, 34, 35,
                                                                       119, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 296, 0, 3, 35, 36,
                                                                       122, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 302, 0, 3, 36, 37,
                                                                       125, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 308, 0, 3, 37, 38,
                                                                       128, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 314, 0, 3, 38, 39,
                                                                       131, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 320, 0, 3, 39, 40,
                                                                       134, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 326, 0, 3, 40, 41,
                                                                       137, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 332, 0, 3, 41, 42,
                                                                       140, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 44, 47,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 47, 50,
                                                                       152, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 50, 53,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 53, 56,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 56, 59,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 59, 62,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 62, 65,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 65, 68,
                                                                       188, 194, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 68, 71,
                                                                       194, 200, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 71, 74,
                                                                       200, 206, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 74, 77,
                                                                       206, 212, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 77, 80,
                                                                       212, 218, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 80, 83,
                                                                       218, 224, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 83, 86,
                                                                       224, 230, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 478, 0, 3, 86, 89,
                                                                       230, 236, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 95, 98,
                                                                       242, 248, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 98,
                                                                       101, 248, 254, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 508, 0, 3, 101,
                                                                       104, 254, 260, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 104,
                                                                       107, 260, 266, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 528, 0, 3, 107,
                                                                       110, 266, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 538, 0, 3, 110,
                                                                       113, 272, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 548, 0, 3, 113,
                                                                       116, 278, 284, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 558, 0, 3, 116,
                                                                       119, 284, 290, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 568, 0, 3, 119,
                                                                       122, 290, 296, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 578, 0, 3, 122,
                                                                       125, 296, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 125,
                                                                       128, 302, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 598, 0, 3, 128,
                                                                       131, 308, 314, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 608, 0, 3, 131,
                                                                       134, 314, 320, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 618, 0, 3, 134,
                                                                       137, 320, 326, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 628, 0, 3, 137,
                                                                       140, 326, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 638, 0, 3, 146,
                                                                       152, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 653, 0, 3, 152,
                                                                       158, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 668, 0, 3, 158,
                                                                       164, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 683, 0, 3, 164,
                                                                       170, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 698, 0, 3, 170,
                                                                       176, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 713, 0, 3, 176,
                                                                       182, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 728, 0, 3, 182,
                                                                       188, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 743, 0, 3, 188,
                                                                       194, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 758, 0, 3, 194,
                                                                       200, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 773, 0, 3, 200,
                                                                       206, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 788, 0, 3, 206,
                                                                       212, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 803, 0, 3, 212,
                                                                       218, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 818, 0, 3, 218,
                                                                       224, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 833, 0, 3, 224,
                                                                       230, 468, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 848, 0, 3, 242,
                                                                       248, 488, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 863, 0, 3, 248,
                                                                       254, 498, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 878, 0, 3, 254,
                                                                       260, 508, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 893, 0, 3, 260,
                                                                       266, 518, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 908, 0, 3, 266,
                                                                       272, 528, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 923, 0, 3, 272,
                                                                       278, 538, 548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 938, 0, 3, 278,
                                                                       284, 548, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 953, 0, 3, 284,
                                                                       290, 558, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 968, 0, 3, 290,
                                                                       296, 568, 578, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 983, 0, 3, 296,
                                                                       302, 578, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 998, 0, 3, 302,
                                                                       308, 588, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1013, 0, 3, 308,
                                                                       314, 598, 608, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1028, 0, 3, 314,
                                                                       320, 608, 618, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 1043, 0, 3, 320,
                                                                       326, 618, 628, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 338,
                                                                       348, 638, 653, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1079, 0, 3, 348,
                                                                       358, 653, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 358,
                                                                       368, 668, 683, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1121, 0, 3, 368,
                                                                       378, 683, 698, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 378,
                                                                       388, 698, 713, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 388,
                                                                       398, 713, 728, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 398,
                                                                       408, 728, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1205, 0, 3, 408,
                                                                       418, 743, 758, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 418,
                                                                       428, 758, 773, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1247, 0, 3, 428,
                                                                       438, 773, 788, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 438,
                                                                       448, 788, 803, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1289, 0, 3, 448,
                                                                       458, 803, 818, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 458,
                                                                       468, 818, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1331, 0, 3, 488,
                                                                       498, 848, 863, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 498,
                                                                       508, 863, 878, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1373, 0, 3, 508,
                                                                       518, 878, 893, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1394, 0, 3, 518,
                                                                       528, 893, 908, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1415, 0, 3, 528,
                                                                       538, 908, 923, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 538,
                                                                       548, 923, 938, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1457, 0, 3, 548,
                                                                       558, 938, 953, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 558,
                                                                       568, 953, 968, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1499, 0, 3, 568,
                                                                       578, 968, 983, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 578,
                                                                       588, 983, 998, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1541, 0, 3, 588,
                                                                       598, 998, 1013, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1562, 0, 3, 598,
                                                                       608, 1013, 1028, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 1583, 0, 3, 608,
                                                                       618, 1028, 1043, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 638,
                                                                       653, 1058, 1079, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1632, 0, 3, 653,
                                                                       668, 1079, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 668,
                                                                       683, 1100, 1121, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 683,
                                                                       698, 1121, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 698,
                                                                       713, 1142, 1163, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 713,
                                                                       728, 1163, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1772, 0, 3, 728,
                                                                       743, 1184, 1205, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 743,
                                                                       758, 1205, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 758,
                                                                       773, 1226, 1247, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 773,
                                                                       788, 1247, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1884, 0, 3, 788,
                                                                       803, 1268, 1289, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1912, 0, 3, 803,
                                                                       818, 1289, 1310, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1940, 0, 3, 848,
                                                                       863, 1331, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1968, 0, 3, 863,
                                                                       878, 1352, 1373, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1996, 0, 3, 878,
                                                                       893, 1373, 1394, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2024, 0, 3, 893,
                                                                       908, 1394, 1415, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 908,
                                                                       923, 1415, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 923,
                                                                       938, 1436, 1457, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2108, 0, 3, 938,
                                                                       953, 1457, 1478, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2136, 0, 3, 953,
                                                                       968, 1478, 1499, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2164, 0, 3, 968,
                                                                       983, 1499, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2192, 0, 3, 983,
                                                                       998, 1520, 1541, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2220, 0, 3, 998,
                                                                       1013, 1541, 1562, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 1013,
                                                                       1028, 1562, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 1058,
                                                                       1079, 1604, 1632, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2312, 0, 3, 1079,
                                                                       1100, 1632, 1660, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2348, 0, 3, 1100,
                                                                       1121, 1660, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2384, 0, 3, 1121,
                                                                       1142, 1688, 1716, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2420, 0, 3, 1142,
                                                                       1163, 1716, 1744, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2456, 0, 3, 1163,
                                                                       1184, 1744, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2492, 0, 3, 1184,
                                                                       1205, 1772, 1800, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2528, 0, 3, 1205,
                                                                       1226, 1800, 1828, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2564, 0, 3, 1226,
                                                                       1247, 1828, 1856, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2600, 0, 3, 1247,
                                                                       1268, 1856, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2636, 0, 3, 1268,
                                                                       1289, 1884, 1912, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2672, 0, 3, 1331,
                                                                       1352, 1940, 1968, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2708, 0, 3, 1352,
                                                                       1373, 1968, 1996, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2744, 0, 3, 1373,
                                                                       1394, 1996, 2024, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2780, 0, 3, 1394,
                                                                       1415, 2024, 2052, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2816, 0, 3, 1415,
                                                                       1436, 2052, 2080, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2852, 0, 3, 1436,
                                                                       1457, 2080, 2108, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2888, 0, 3, 1457,
                                                                       1478, 2108, 2136, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2924, 0, 3, 1478,
                                                                       1499, 2136, 2164, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2960, 0, 3, 1499,
                                                                       1520, 2164, 2192, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 2996, 0, 3, 1520,
                                                                       1541, 2192, 2220, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 3032, 0, 3, 1541,
                                                                       1562, 2220, 2248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3068, 0, 3, 1604,
                                                                       1632, 2276, 2312, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3113, 0, 3, 1632,
                                                                       1660, 2312, 2348, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3158, 0, 3, 1660,
                                                                       1688, 2348, 2384, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3203, 0, 3, 1688,
                                                                       1716, 2384, 2420, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3248, 0, 3, 1716,
                                                                       1744, 2420, 2456, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3293, 0, 3, 1744,
                                                                       1772, 2456, 2492, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3338, 0, 3, 1772,
                                                                       1800, 2492, 2528, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3383, 0, 3, 1800,
                                                                       1828, 2528, 2564, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3428, 0, 3, 1828,
                                                                       1856, 2564, 2600, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3473, 0, 3, 1856,
                                                                       1884, 2600, 2636, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3518, 0, 3, 1940,
                                                                       1968, 2672, 2708, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3563, 0, 3, 1968,
                                                                       1996, 2708, 2744, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3608, 0, 3, 1996,
                                                                       2024, 2744, 2780, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3653, 0, 3, 2024,
                                                                       2052, 2780, 2816, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3698, 0, 3, 2052,
                                                                       2080, 2816, 2852, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3743, 0, 3, 2080,
                                                                       2108, 2852, 2888, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3788, 0, 3, 2108,
                                                                       2136, 2888, 2924, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3833, 0, 3, 2136,
                                                                       2164, 2924, 2960, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3878, 0, 3, 2164,
                                                                       2192, 2960, 2996, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 3923, 0, 3, 2192,
                                                                       2220, 2996, 3032, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 3968, 0, 3, 2276,
                                                                       2312, 3068, 3113, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4023, 0, 3, 2312,
                                                                       2348, 3113, 3158, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4078, 0, 3, 2348,
                                                                       2384, 3158, 3203, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4133, 0, 3, 2384,
                                                                       2420, 3203, 3248, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4188, 0, 3, 2420,
                                                                       2456, 3248, 3293, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4243, 0, 3, 2456,
                                                                       2492, 3293, 3338, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 2492,
                                                                       2528, 3338, 3383, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4353, 0, 3, 2528,
                                                                       2564, 3383, 3428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4408, 0, 3, 2564,
                                                                       2600, 3428, 3473, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4463, 0, 3, 2672,
                                                                       2708, 3518, 3563, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4518, 0, 3, 2708,
                                                                       2744, 3563, 3608, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4573, 0, 3, 2744,
                                                                       2780, 3608, 3653, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 2780,
                                                                       2816, 3653, 3698, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4683, 0, 3, 2816,
                                                                       2852, 3698, 3743, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4738, 0, 3, 2852,
                                                                       2888, 3743, 3788, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4793, 0, 3, 2888,
                                                                       2924, 3788, 3833, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4848, 0, 3, 2924,
                                                                       2960, 3833, 3878, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 4903, 0, 3, 2960,
                                                                       2996, 3878, 3923, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 4958, 0, 3, 3068,
                                                                       3113, 3968, 4023, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5024, 0, 3, 3113,
                                                                       3158, 4023, 4078, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5090, 0, 3, 3158,
                                                                       3203, 4078, 4133, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5156, 0, 3, 3203,
                                                                       3248, 4133, 4188, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5222, 0, 3, 3248,
                                                                       3293, 4188, 4243, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5288, 0, 3, 3293,
                                                                       3338, 4243, 4298, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5354, 0, 3, 3338,
                                                                       3383, 4298, 4353, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5420, 0, 3, 3383,
                                                                       3428, 4353, 4408, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5486, 0, 3, 3518,
                                                                       3563, 4463, 4518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5552, 0, 3, 3563,
                                                                       3608, 4518, 4573, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5618, 0, 3, 3608,
                                                                       3653, 4573, 4628, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5684, 0, 3, 3653,
                                                                       3698, 4628, 4683, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5750, 0, 3, 3698,
                                                                       3743, 4683, 4738, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5816, 0, 3, 3743,
                                                                       3788, 4738, 4793, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5882, 0, 3, 3788,
                                                                       3833, 4793, 4848, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 5948, 0, 3, 3833,
                                                                       3878, 4848, 4903, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6014, 0, 3, 3968,
                                                                       4023, 4958, 5024, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6092, 0, 3, 4023,
                                                                       4078, 5024, 5090, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6170, 0, 3, 4078,
                                                                       4133, 5090, 5156, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 4133,
                                                                       4188, 5156, 5222, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6326, 0, 3, 4188,
                                                                       4243, 5222, 5288, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6404, 0, 3, 4243,
                                                                       4298, 5288, 5354, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6482, 0, 3, 4298,
                                                                       4353, 5354, 5420, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6560, 0, 3, 4463,
                                                                       4518, 5486, 5552, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6638, 0, 3, 4518,
                                                                       4573, 5552, 5618, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 4573,
                                                                       4628, 5618, 5684, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6794, 0, 3, 4628,
                                                                       4683, 5684, 5750, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6872, 0, 3, 4683,
                                                                       4738, 5750, 5816, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 6950, 0, 3, 4738,
                                                                       4793, 5816, 5882, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 7028, 0, 3, 4793,
                                                                       4848, 5882, 5948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7106, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7109, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7112, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7115, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7118, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7121, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7124, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7127, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7130, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7133, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7136, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7139, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7142, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7145, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7148, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7151, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7154, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7157, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7160, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7163, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7166, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7169, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7172, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7175, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7178, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7181, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7184, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7187, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7190, 3, 36,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7193, 3, 37,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7196, 3, 38,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7199, 3, 39,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7202, 3, 40,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7205, 3, 41,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7208, 3, 42,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 7211, 3, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7214, 3, 9, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7223, 3, 10, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7232, 3, 11, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7241, 3, 12, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7250, 3, 13, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7259, 3, 14, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7268, 3, 15, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7277, 3, 16, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7286, 3, 17, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7295, 3, 18, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7304, 3, 19, 80,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7313, 3, 20, 83,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7322, 3, 21, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7331, 3, 22, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7340, 3, 23, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7349, 3, 28, 101,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7358, 3, 29, 104,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7367, 3, 30, 107,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7376, 3, 31, 110,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7385, 3, 32, 113,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7394, 3, 33, 116,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7403, 3, 34, 119,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7412, 3, 35, 122,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7421, 3, 36, 125,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7430, 3, 37, 128,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7439, 3, 38, 131,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7448, 3, 39, 134,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7457, 3, 40, 137,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7466, 3, 41, 140,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 7475, 3, 42, 143,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7484, 3, 44, 146,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7502, 3, 47, 152,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7520, 3, 50, 158,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7538, 3, 53, 164,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7556, 3, 56, 170,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7574, 3, 59, 176,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7592, 3, 62, 182,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7610, 3, 65, 188,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7628, 3, 68, 194,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7646, 3, 71, 200,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7664, 3, 74, 206,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7682, 3, 77, 212,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7700, 3, 80, 218,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7718, 3, 83, 224,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7736, 3, 86, 230,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7754, 3, 89, 236,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7772, 3, 95, 242,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7790, 3, 98, 248,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7808, 3, 101, 254,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7826, 3, 104, 260,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7844, 3, 107, 266,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7862, 3, 110, 272,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7880, 3, 113, 278,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7898, 3, 116, 284,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7916, 3, 119, 290,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7934, 3, 122, 296,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7952, 3, 125, 302,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7970, 3, 128, 308,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 7988, 3, 131, 314,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 8006, 3, 134, 320,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 8024, 3, 137, 326,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 8042, 3, 140, 332,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8060, 3, 146, 338,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8090, 3, 152, 348,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8120, 3, 158, 358,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8150, 3, 164, 368,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8180, 3, 170, 378,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8210, 3, 176, 388,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8240, 3, 182, 398,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8270, 3, 188, 408,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8300, 3, 194, 418,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8330, 3, 200, 428,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8360, 3, 206, 438,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8390, 3, 212, 448,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8420, 3, 218, 458,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8450, 3, 224, 468,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8480, 3, 230, 478,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8510, 3, 242, 488,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8540, 3, 248, 498,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8570, 3, 254, 508,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8600, 3, 260, 518,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8630, 3, 266, 528,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8660, 3, 272, 538,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8690, 3, 278, 548,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8720, 3, 284, 558,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8750, 3, 290, 568,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8780, 3, 296, 578,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8810, 3, 302, 588,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8840, 3, 308, 598,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8870, 3, 314, 608,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8900, 3, 320, 618,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 8930, 3, 326, 628,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 8960, 3, 338, 638,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9005, 3, 348, 653,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9050, 3, 358, 668,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9095, 3, 368, 683,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9140, 3, 378, 698,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9185, 3, 388, 713,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9230, 3, 398, 728,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9275, 3, 408, 743,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9320, 3, 418, 758,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9365, 3, 428, 773,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9410, 3, 438, 788,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9455, 3, 448, 803,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9500, 3, 458, 818,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9545, 3, 468, 833,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9590, 3, 488, 848,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9635, 3, 498, 863,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9680, 3, 508, 878,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9725, 3, 518, 893,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9770, 3, 528, 908,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9815, 3, 538, 923,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9860, 3, 548, 938,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9905, 3, 558, 953,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9950, 3, 568, 968,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 9995, 3, 578, 983,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10040, 3, 588,
                                                                       998, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10085, 3, 598,
                                                                       1013, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10130, 3, 608,
                                                                       1028, ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 10175, 3, 618,
                                                                       1043, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10220, 3, 638,
                                                                       1058, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10283, 3, 653,
                                                                       1079, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10346, 3, 668,
                                                                       1100, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10409, 3, 683,
                                                                       1121, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10472, 3, 698,
                                                                       1142, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10535, 3, 713,
                                                                       1163, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10598, 3, 728,
                                                                       1184, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10661, 3, 743,
                                                                       1205, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10724, 3, 758,
                                                                       1226, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10787, 3, 773,
                                                                       1247, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10850, 3, 788,
                                                                       1268, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10913, 3, 803,
                                                                       1289, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 10976, 3, 818,
                                                                       1310, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11039, 3, 848,
                                                                       1331, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11102, 3, 863,
                                                                       1352, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11165, 3, 878,
                                                                       1373, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11228, 3, 893,
                                                                       1394, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11291, 3, 908,
                                                                       1415, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11354, 3, 923,
                                                                       1436, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11417, 3, 938,
                                                                       1457, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11480, 3, 953,
                                                                       1478, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11543, 3, 968,
                                                                       1499, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11606, 3, 983,
                                                                       1520, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11669, 3, 998,
                                                                       1541, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11732, 3, 1013,
                                                                       1562, ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 11795, 3, 1028,
                                                                       1583, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11858, 3, 1058,
                                                                       1604, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 11942, 3, 1079,
                                                                       1632, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12026, 3, 1100,
                                                                       1660, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12110, 3, 1121,
                                                                       1688, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12194, 3, 1142,
                                                                       1716, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12278, 3, 1163,
                                                                       1744, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12362, 3, 1184,
                                                                       1772, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12446, 3, 1205,
                                                                       1800, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12530, 3, 1226,
                                                                       1828, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12614, 3, 1247,
                                                                       1856, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12698, 3, 1268,
                                                                       1884, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12782, 3, 1289,
                                                                       1912, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12866, 3, 1331,
                                                                       1940, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 12950, 3, 1352,
                                                                       1968, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13034, 3, 1373,
                                                                       1996, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13118, 3, 1394,
                                                                       2024, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13202, 3, 1415,
                                                                       2052, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13286, 3, 1436,
                                                                       2080, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13370, 3, 1457,
                                                                       2108, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13454, 3, 1478,
                                                                       2136, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13538, 3, 1499,
                                                                       2164, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13622, 3, 1520,
                                                                       2192, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13706, 3, 1541,
                                                                       2220, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 13790, 3, 1562,
                                                                       2248, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13874, 3, 1604,
                                                                       2276, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 13982, 3, 1632,
                                                                       2312, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14090, 3, 1660,
                                                                       2348, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14198, 3, 1688,
                                                                       2384, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14306, 3, 1716,
                                                                       2420, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14414, 3, 1744,
                                                                       2456, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14522, 3, 1772,
                                                                       2492, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14630, 3, 1800,
                                                                       2528, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14738, 3, 1828,
                                                                       2564, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14846, 3, 1856,
                                                                       2600, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 14954, 3, 1884,
                                                                       2636, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15062, 3, 1940,
                                                                       2672, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15170, 3, 1968,
                                                                       2708, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15278, 3, 1996,
                                                                       2744, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15386, 3, 2024,
                                                                       2780, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15494, 3, 2052,
                                                                       2816, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15602, 3, 2080,
                                                                       2852, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15710, 3, 2108,
                                                                       2888, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15818, 3, 2136,
                                                                       2924, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 15926, 3, 2164,
                                                                       2960, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16034, 3, 2192,
                                                                       2996, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 16142, 3, 2220,
                                                                       3032, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16250, 3, 2276,
                                                                       3068, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16385, 3, 2312,
                                                                       3113, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16520, 3, 2348,
                                                                       3158, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16655, 3, 2384,
                                                                       3203, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16790, 3, 2420,
                                                                       3248, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 16925, 3, 2456,
                                                                       3293, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17060, 3, 2492,
                                                                       3338, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17195, 3, 2528,
                                                                       3383, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17330, 3, 2564,
                                                                       3428, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17465, 3, 2600,
                                                                       3473, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17600, 3, 2672,
                                                                       3518, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17735, 3, 2708,
                                                                       3563, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 17870, 3, 2744,
                                                                       3608, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18005, 3, 2780,
                                                                       3653, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18140, 3, 2816,
                                                                       3698, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18275, 3, 2852,
                                                                       3743, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18410, 3, 2888,
                                                                       3788, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18545, 3, 2924,
                                                                       3833, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18680, 3, 2960,
                                                                       3878, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 18815, 3, 2996,
                                                                       3923, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 18950, 3, 3068,
                                                                       3968, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19115, 3, 3113,
                                                                       4023, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19280, 3, 3158,
                                                                       4078, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19445, 3, 3203,
                                                                       4133, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19610, 3, 3248,
                                                                       4188, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19775, 3, 3293,
                                                                       4243, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 19940, 3, 3338,
                                                                       4298, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 20105, 3, 3383,
                                                                       4353, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 20270, 3, 3428,
                                                                       4408, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 20435, 3, 3518,
                                                                       4463, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 20600, 3, 3563,
                                                                       4518, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 20765, 3, 3608,
                                                                       4573, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 20930, 3, 3653,
                                                                       4628, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21095, 3, 3698,
                                                                       4683, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21260, 3, 3743,
                                                                       4738, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21425, 3, 3788,
                                                                       4793, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21590, 3, 3833,
                                                                       4848, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 21755, 3, 3878,
                                                                       4903, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 21920, 3, 3968,
                                                                       4958, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 22118, 3, 4023,
                                                                       5024, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 22316, 3, 4078,
                                                                       5090, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 22514, 3, 4133,
                                                                       5156, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 22712, 3, 4188,
                                                                       5222, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 22910, 3, 4243,
                                                                       5288, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 23108, 3, 4298,
                                                                       5354, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 23306, 3, 4353,
                                                                       5420, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 23504, 3, 4463,
                                                                       5486, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 23702, 3, 4518,
                                                                       5552, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 23900, 3, 4573,
                                                                       5618, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24098, 3, 4628,
                                                                       5684, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24296, 3, 4683,
                                                                       5750, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24494, 3, 4738,
                                                                       5816, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24692, 3, 4793,
                                                                       5882, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 24890, 3, 4848,
                                                                       5948, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 25088, 3, 4958,
                                                                       6014, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 25322, 3, 5024,
                                                                       6092, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 25556, 3, 5090,
                                                                       6170, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 25790, 3, 5156,
                                                                       6248, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 26024, 3, 5222,
                                                                       6326, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 26258, 3, 5288,
                                                                       6404, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 26492, 3, 5354,
                                                                       6482, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 26726, 3, 5486,
                                                                       6560, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 26960, 3, 5552,
                                                                       6638, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 27194, 3, 5618,
                                                                       6716, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 27428, 3, 5684,
                                                                       6794, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 27662, 3, 5750,
                                                                       6872, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 27896, 3, 5816,
                                                                       6950, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 28130, 3, 5882,
                                                                       7028, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28364, 3, 7, 8,
                                                                       7112, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28370, 3, 8, 9,
                                                                       7115, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28376, 3, 9, 10,
                                                                       7118, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28382, 3, 10, 11,
                                                                       7121, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28388, 3, 11, 12,
                                                                       7124, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28394, 3, 12, 13,
                                                                       7127, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28400, 3, 13, 14,
                                                                       7130, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28406, 3, 14, 15,
                                                                       7133, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28412, 3, 15, 16,
                                                                       7136, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28418, 3, 16, 17,
                                                                       7139, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28424, 3, 17, 18,
                                                                       7142, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28430, 3, 18, 19,
                                                                       7145, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28436, 3, 19, 20,
                                                                       7148, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28442, 3, 20, 21,
                                                                       7151, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28448, 3, 21, 22,
                                                                       7154, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28454, 3, 22, 23,
                                                                       7157, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28460, 3, 26, 27,
                                                                       7166, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28466, 3, 27, 28,
                                                                       7169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28472, 3, 28, 29,
                                                                       7172, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28478, 3, 29, 30,
                                                                       7175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28484, 3, 30, 31,
                                                                       7178, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28490, 3, 31, 32,
                                                                       7181, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28496, 3, 32, 33,
                                                                       7184, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28502, 3, 33, 34,
                                                                       7187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28508, 3, 34, 35,
                                                                       7190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28514, 3, 35, 36,
                                                                       7193, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28520, 3, 36, 37,
                                                                       7196, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28526, 3, 37, 38,
                                                                       7199, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28532, 3, 38, 39,
                                                                       7202, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28538, 3, 39, 40,
                                                                       7205, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28544, 3, 40, 41,
                                                                       7208, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 28550, 3, 41, 42,
                                                                       7211, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28556, 0, 3,
                                                                       28364, 7112, 28370, 7214,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28574, 0, 3,
                                                                       28370, 7115, 28376, 7223,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28592, 0, 3,
                                                                       28376, 7118, 28382, 7232,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28610, 0, 3,
                                                                       28382, 7121, 28388, 7241,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28628, 0, 3,
                                                                       28388, 7124, 28394, 7250,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28646, 0, 3,
                                                                       28394, 7127, 28400, 7259,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28664, 0, 3,
                                                                       28400, 7130, 28406, 7268,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28682, 0, 3,
                                                                       28406, 7133, 28412, 7277,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28700, 0, 3,
                                                                       28412, 7136, 28418, 7286,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28718, 0, 3,
                                                                       28418, 7139, 28424, 7295,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28736, 0, 3,
                                                                       28424, 7142, 28430, 7304,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28754, 0, 3,
                                                                       28430, 7145, 28436, 7313,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28772, 0, 3,
                                                                       28436, 7148, 28442, 7322,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28790, 0, 3,
                                                                       28442, 7151, 28448, 7331,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28808, 0, 3,
                                                                       28448, 7154, 28454, 7340,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28826, 0, 3,
                                                                       28460, 7166, 28466, 7349,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28844, 0, 3,
                                                                       28466, 7169, 28472, 7358,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28862, 0, 3,
                                                                       28472, 7172, 28478, 7367,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28880, 0, 3,
                                                                       28478, 7175, 28484, 7376,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28898, 0, 3,
                                                                       28484, 7178, 28490, 7385,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28916, 0, 3,
                                                                       28490, 7181, 28496, 7394,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28934, 0, 3,
                                                                       28496, 7184, 28502, 7403,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28952, 0, 3,
                                                                       28502, 7187, 28508, 7412,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28970, 0, 3,
                                                                       28508, 7190, 28514, 7421,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 28988, 0, 3,
                                                                       28514, 7193, 28520, 7430,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 29006, 0, 3,
                                                                       28520, 7196, 28526, 7439,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 29024, 0, 3,
                                                                       28526, 7199, 28532, 7448,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 29042, 0, 3,
                                                                       28532, 7202, 28538, 7457,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 29060, 0, 3,
                                                                       28538, 7205, 28544, 7466,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 29078, 0, 3,
                                                                       28544, 7208, 28550, 7475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29096, 0, 3,
                                                                       28556, 7214, 28574, 146,
                                                                       152, 7520, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29132, 0, 3,
                                                                       28574, 7223, 28592, 152,
                                                                       158, 7538, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29168, 0, 3,
                                                                       28592, 7232, 28610, 158,
                                                                       164, 7556, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29204, 0, 3,
                                                                       28610, 7241, 28628, 164,
                                                                       170, 7574, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29240, 0, 3,
                                                                       28628, 7250, 28646, 170,
                                                                       176, 7592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29276, 0, 3,
                                                                       28646, 7259, 28664, 176,
                                                                       182, 7610, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29312, 0, 3,
                                                                       28664, 7268, 28682, 182,
                                                                       188, 7628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29348, 0, 3,
                                                                       28682, 7277, 28700, 188,
                                                                       194, 7646, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29384, 0, 3,
                                                                       28700, 7286, 28718, 194,
                                                                       200, 7664, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29420, 0, 3,
                                                                       28718, 7295, 28736, 200,
                                                                       206, 7682, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29456, 0, 3,
                                                                       28736, 7304, 28754, 206,
                                                                       212, 7700, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29492, 0, 3,
                                                                       28754, 7313, 28772, 212,
                                                                       218, 7718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29528, 0, 3,
                                                                       28772, 7322, 28790, 218,
                                                                       224, 7736, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29564, 0, 3,
                                                                       28790, 7331, 28808, 224,
                                                                       230, 7754, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29600, 0, 3,
                                                                       28826, 7349, 28844, 242,
                                                                       248, 7808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29636, 0, 3,
                                                                       28844, 7358, 28862, 248,
                                                                       254, 7826, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29672, 0, 3,
                                                                       28862, 7367, 28880, 254,
                                                                       260, 7844, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29708, 0, 3,
                                                                       28880, 7376, 28898, 260,
                                                                       266, 7862, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29744, 0, 3,
                                                                       28898, 7385, 28916, 266,
                                                                       272, 7880, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29780, 0, 3,
                                                                       28916, 7394, 28934, 272,
                                                                       278, 7898, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29816, 0, 3,
                                                                       28934, 7403, 28952, 278,
                                                                       284, 7916, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29852, 0, 3,
                                                                       28952, 7412, 28970, 284,
                                                                       290, 7934, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29888, 0, 3,
                                                                       28970, 7421, 28988, 290,
                                                                       296, 7952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29924, 0, 3,
                                                                       28988, 7430, 29006, 296,
                                                                       302, 7970, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29960, 0, 3,
                                                                       29006, 7439, 29024, 302,
                                                                       308, 7988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 29996, 0, 3,
                                                                       29024, 7448, 29042, 308,
                                                                       314, 8006, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30032, 0, 3,
                                                                       29042, 7457, 29060, 314,
                                                                       320, 8024, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 30068, 0, 3,
                                                                       29060, 7466, 29078, 320,
                                                                       326, 8042, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30104, 0, 3,
                                                                       29096, 7520, 29132, 338,
                                                                       348, 8120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30164, 0, 3,
                                                                       29132, 7538, 29168, 348,
                                                                       358, 8150, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30224, 0, 3,
                                                                       29168, 7556, 29204, 358,
                                                                       368, 8180, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30284, 0, 3,
                                                                       29204, 7574, 29240, 368,
                                                                       378, 8210, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30344, 0, 3,
                                                                       29240, 7592, 29276, 378,
                                                                       388, 8240, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30404, 0, 3,
                                                                       29276, 7610, 29312, 388,
                                                                       398, 8270, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30464, 0, 3,
                                                                       29312, 7628, 29348, 398,
                                                                       408, 8300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30524, 0, 3,
                                                                       29348, 7646, 29384, 408,
                                                                       418, 8330, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30584, 0, 3,
                                                                       29384, 7664, 29420, 418,
                                                                       428, 8360, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30644, 0, 3,
                                                                       29420, 7682, 29456, 428,
                                                                       438, 8390, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30704, 0, 3,
                                                                       29456, 7700, 29492, 438,
                                                                       448, 8420, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30764, 0, 3,
                                                                       29492, 7718, 29528, 448,
                                                                       458, 8450, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30824, 0, 3,
                                                                       29528, 7736, 29564, 458,
                                                                       468, 8480, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30884, 0, 3,
                                                                       29600, 7808, 29636, 488,
                                                                       498, 8570, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 30944, 0, 3,
                                                                       29636, 7826, 29672, 498,
                                                                       508, 8600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31004, 0, 3,
                                                                       29672, 7844, 29708, 508,
                                                                       518, 8630, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31064, 0, 3,
                                                                       29708, 7862, 29744, 518,
                                                                       528, 8660, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31124, 0, 3,
                                                                       29744, 7880, 29780, 528,
                                                                       538, 8690, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31184, 0, 3,
                                                                       29780, 7898, 29816, 538,
                                                                       548, 8720, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31244, 0, 3,
                                                                       29816, 7916, 29852, 548,
                                                                       558, 8750, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31304, 0, 3,
                                                                       29852, 7934, 29888, 558,
                                                                       568, 8780, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31364, 0, 3,
                                                                       29888, 7952, 29924, 568,
                                                                       578, 8810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31424, 0, 3,
                                                                       29924, 7970, 29960, 578,
                                                                       588, 8840, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31484, 0, 3,
                                                                       29960, 7988, 29996, 588,
                                                                       598, 8870, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31544, 0, 3,
                                                                       29996, 8006, 30032, 598,
                                                                       608, 8900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 31604, 0, 3,
                                                                       30032, 8024, 30068, 608,
                                                                       618, 8930, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 31664, 0, 3,
                                                                       30104, 8120, 30164, 638,
                                                                       653, 9050, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 31754, 0, 3,
                                                                       30164, 8150, 30224, 653,
                                                                       668, 9095, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 31844, 0, 3,
                                                                       30224, 8180, 30284, 668,
                                                                       683, 9140, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 31934, 0, 3,
                                                                       30284, 8210, 30344, 683,
                                                                       698, 9185, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32024, 0, 3,
                                                                       30344, 8240, 30404, 698,
                                                                       713, 9230, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32114, 0, 3,
                                                                       30404, 8270, 30464, 713,
                                                                       728, 9275, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32204, 0, 3,
                                                                       30464, 8300, 30524, 728,
                                                                       743, 9320, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32294, 0, 3,
                                                                       30524, 8330, 30584, 743,
                                                                       758, 9365, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32384, 0, 3,
                                                                       30584, 8360, 30644, 758,
                                                                       773, 9410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32474, 0, 3,
                                                                       30644, 8390, 30704, 773,
                                                                       788, 9455, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32564, 0, 3,
                                                                       30704, 8420, 30764, 788,
                                                                       803, 9500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32654, 0, 3,
                                                                       30764, 8450, 30824, 803,
                                                                       818, 9545, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32744, 0, 3,
                                                                       30884, 8570, 30944, 848,
                                                                       863, 9680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32834, 0, 3,
                                                                       30944, 8600, 31004, 863,
                                                                       878, 9725, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 32924, 0, 3,
                                                                       31004, 8630, 31064, 878,
                                                                       893, 9770, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33014, 0, 3,
                                                                       31064, 8660, 31124, 893,
                                                                       908, 9815, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33104, 0, 3,
                                                                       31124, 8690, 31184, 908,
                                                                       923, 9860, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33194, 0, 3,
                                                                       31184, 8720, 31244, 923,
                                                                       938, 9905, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33284, 0, 3,
                                                                       31244, 8750, 31304, 938,
                                                                       953, 9950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33374, 0, 3,
                                                                       31304, 8780, 31364, 953,
                                                                       968, 9995, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33464, 0, 3,
                                                                       31364, 8810, 31424, 968,
                                                                       983, 10040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33554, 0, 3,
                                                                       31424, 8840, 31484, 983,
                                                                       998, 10085, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33644, 0, 3,
                                                                       31484, 8870, 31544, 998,
                                                                       1013, 10130, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 33734, 0, 3,
                                                                       31544, 8900, 31604, 1013,
                                                                       1028, 10175, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 33824, 0, 3,
                                                                       31664, 9050, 31754, 1058,
                                                                       1079, 10346, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 33950, 0, 3,
                                                                       31754, 9095, 31844, 1079,
                                                                       1100, 10409, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 34076, 0, 3,
                                                                       31844, 9140, 31934, 1100,
                                                                       1121, 10472, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 34202, 0, 3,
                                                                       31934, 9185, 32024, 1121,
                                                                       1142, 10535, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 34328, 0, 3,
                                                                       32024, 9230, 32114, 1142,
                                                                       1163, 10598, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 34454, 0, 3,
                                                                       32114, 9275, 32204, 1163,
                                                                       1184, 10661, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 34580, 0, 3,
                                                                       32204, 9320, 32294, 1184,
                                                                       1205, 10724, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 34706, 0, 3,
                                                                       32294, 9365, 32384, 1205,
                                                                       1226, 10787, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 34832, 0, 3,
                                                                       32384, 9410, 32474, 1226,
                                                                       1247, 10850, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 34958, 0, 3,
                                                                       32474, 9455, 32564, 1247,
                                                                       1268, 10913, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 35084, 0, 3,
                                                                       32564, 9500, 32654, 1268,
                                                                       1289, 10976, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 35210, 0, 3,
                                                                       32744, 9680, 32834, 1331,
                                                                       1352, 11165, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 35336, 0, 3,
                                                                       32834, 9725, 32924, 1352,
                                                                       1373, 11228, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 35462, 0, 3,
                                                                       32924, 9770, 33014, 1373,
                                                                       1394, 11291, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 35588, 0, 3,
                                                                       33014, 9815, 33104, 1394,
                                                                       1415, 11354, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 35714, 0, 3,
                                                                       33104, 9860, 33194, 1415,
                                                                       1436, 11417, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 35840, 0, 3,
                                                                       33194, 9905, 33284, 1436,
                                                                       1457, 11480, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 35966, 0, 3,
                                                                       33284, 9950, 33374, 1457,
                                                                       1478, 11543, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36092, 0, 3,
                                                                       33374, 9995, 33464, 1478,
                                                                       1499, 11606, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36218, 0, 3,
                                                                       33464, 10040, 33554, 1499,
                                                                       1520, 11669, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36344, 0, 3,
                                                                       33554, 10085, 33644, 1520,
                                                                       1541, 11732, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 36470, 0, 3,
                                                                       33644, 10130, 33734, 1541,
                                                                       1562, 11795, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 36596, 0, 3,
                                                                       33824, 10346, 33950, 1604,
                                                                       1632, 12026, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 36764, 0, 3,
                                                                       33950, 10409, 34076, 1632,
                                                                       1660, 12110, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 36932, 0, 3,
                                                                       34076, 10472, 34202, 1660,
                                                                       1688, 12194, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 37100, 0, 3,
                                                                       34202, 10535, 34328, 1688,
                                                                       1716, 12278, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 37268, 0, 3,
                                                                       34328, 10598, 34454, 1716,
                                                                       1744, 12362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 37436, 0, 3,
                                                                       34454, 10661, 34580, 1744,
                                                                       1772, 12446, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 37604, 0, 3,
                                                                       34580, 10724, 34706, 1772,
                                                                       1800, 12530, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 37772, 0, 3,
                                                                       34706, 10787, 34832, 1800,
                                                                       1828, 12614, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 37940, 0, 3,
                                                                       34832, 10850, 34958, 1828,
                                                                       1856, 12698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 38108, 0, 3,
                                                                       34958, 10913, 35084, 1856,
                                                                       1884, 12782, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 38276, 0, 3,
                                                                       35210, 11165, 35336, 1940,
                                                                       1968, 13034, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 38444, 0, 3,
                                                                       35336, 11228, 35462, 1968,
                                                                       1996, 13118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 38612, 0, 3,
                                                                       35462, 11291, 35588, 1996,
                                                                       2024, 13202, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 38780, 0, 3,
                                                                       35588, 11354, 35714, 2024,
                                                                       2052, 13286, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 38948, 0, 3,
                                                                       35714, 11417, 35840, 2052,
                                                                       2080, 13370, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 39116, 0, 3,
                                                                       35840, 11480, 35966, 2080,
                                                                       2108, 13454, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 39284, 0, 3,
                                                                       35966, 11543, 36092, 2108,
                                                                       2136, 13538, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 39452, 0, 3,
                                                                       36092, 11606, 36218, 2136,
                                                                       2164, 13622, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 39620, 0, 3,
                                                                       36218, 11669, 36344, 2164,
                                                                       2192, 13706, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 39788, 0, 3,
                                                                       36344, 11732, 36470, 2192,
                                                                       2220, 13790, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 39956, 0, 3,
                                                                       36596, 12026, 36764, 2276,
                                                                       2312, 14090, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 40172, 0, 3,
                                                                       36764, 12110, 36932, 2312,
                                                                       2348, 14198, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 40388, 0, 3,
                                                                       36932, 12194, 37100, 2348,
                                                                       2384, 14306, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 40604, 0, 3,
                                                                       37100, 12278, 37268, 2384,
                                                                       2420, 14414, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 40820, 0, 3,
                                                                       37268, 12362, 37436, 2420,
                                                                       2456, 14522, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 41036, 0, 3,
                                                                       37436, 12446, 37604, 2456,
                                                                       2492, 14630, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 41252, 0, 3,
                                                                       37604, 12530, 37772, 2492,
                                                                       2528, 14738, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 41468, 0, 3,
                                                                       37772, 12614, 37940, 2528,
                                                                       2564, 14846, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 41684, 0, 3,
                                                                       37940, 12698, 38108, 2564,
                                                                       2600, 14954, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 41900, 0, 3,
                                                                       38276, 13034, 38444, 2672,
                                                                       2708, 15278, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 42116, 0, 3,
                                                                       38444, 13118, 38612, 2708,
                                                                       2744, 15386, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 42332, 0, 3,
                                                                       38612, 13202, 38780, 2744,
                                                                       2780, 15494, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 42548, 0, 3,
                                                                       38780, 13286, 38948, 2780,
                                                                       2816, 15602, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 42764, 0, 3,
                                                                       38948, 13370, 39116, 2816,
                                                                       2852, 15710, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 42980, 0, 3,
                                                                       39116, 13454, 39284, 2852,
                                                                       2888, 15818, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 43196, 0, 3,
                                                                       39284, 13538, 39452, 2888,
                                                                       2924, 15926, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 43412, 0, 3,
                                                                       39452, 13622, 39620, 2924,
                                                                       2960, 16034, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 43628, 0, 3,
                                                                       39620, 13706, 39788, 2960,
                                                                       2996, 16142, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 43844, 0, 3,
                                                                       39956, 14090, 40172, 3068,
                                                                       3113, 16520, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 44114, 0, 3,
                                                                       40172, 14198, 40388, 3113,
                                                                       3158, 16655, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 44384, 0, 3,
                                                                       40388, 14306, 40604, 3158,
                                                                       3203, 16790, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 44654, 0, 3,
                                                                       40604, 14414, 40820, 3203,
                                                                       3248, 16925, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 44924, 0, 3,
                                                                       40820, 14522, 41036, 3248,
                                                                       3293, 17060, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 45194, 0, 3,
                                                                       41036, 14630, 41252, 3293,
                                                                       3338, 17195, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 45464, 0, 3,
                                                                       41252, 14738, 41468, 3338,
                                                                       3383, 17330, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 45734, 0, 3,
                                                                       41468, 14846, 41684, 3383,
                                                                       3428, 17465, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 46004, 0, 3,
                                                                       41900, 15278, 42116, 3518,
                                                                       3563, 17870, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 46274, 0, 3,
                                                                       42116, 15386, 42332, 3563,
                                                                       3608, 18005, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 46544, 0, 3,
                                                                       42332, 15494, 42548, 3608,
                                                                       3653, 18140, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 46814, 0, 3,
                                                                       42548, 15602, 42764, 3653,
                                                                       3698, 18275, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 47084, 0, 3,
                                                                       42764, 15710, 42980, 3698,
                                                                       3743, 18410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 47354, 0, 3,
                                                                       42980, 15818, 43196, 3743,
                                                                       3788, 18545, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 47624, 0, 3,
                                                                       43196, 15926, 43412, 3788,
                                                                       3833, 18680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 47894, 0, 3,
                                                                       43412, 16034, 43628, 3833,
                                                                       3878, 18815, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 48164, 0, 3,
                                                                       43844, 16520, 44114, 3968,
                                                                       4023, 19280, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 48494, 0, 3,
                                                                       44114, 16655, 44384, 4023,
                                                                       4078, 19445, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 48824, 0, 3,
                                                                       44384, 16790, 44654, 4078,
                                                                       4133, 19610, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 49154, 0, 3,
                                                                       44654, 16925, 44924, 4133,
                                                                       4188, 19775, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 49484, 0, 3,
                                                                       44924, 17060, 45194, 4188,
                                                                       4243, 19940, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 49814, 0, 3,
                                                                       45194, 17195, 45464, 4243,
                                                                       4298, 20105, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 50144, 0, 3,
                                                                       45464, 17330, 45734, 4298,
                                                                       4353, 20270, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 50474, 0, 3,
                                                                       46004, 17870, 46274, 4463,
                                                                       4518, 20765, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 50804, 0, 3,
                                                                       46274, 18005, 46544, 4518,
                                                                       4573, 20930, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 51134, 0, 3,
                                                                       46544, 18140, 46814, 4573,
                                                                       4628, 21095, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 51464, 0, 3,
                                                                       46814, 18275, 47084, 4628,
                                                                       4683, 21260, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 51794, 0, 3,
                                                                       47084, 18410, 47354, 4683,
                                                                       4738, 21425, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 52124, 0, 3,
                                                                       47354, 18545, 47624, 4738,
                                                                       4793, 21590, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 52454, 0, 3,
                                                                       47624, 18680, 47894, 4793,
                                                                       4848, 21755, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 52784, 0, 3,
                                                                       48164, 19280, 48494, 4958,
                                                                       5024, 22316, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 53180, 0, 3,
                                                                       48494, 19445, 48824, 5024,
                                                                       5090, 22514, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 53576, 0, 3,
                                                                       48824, 19610, 49154, 5090,
                                                                       5156, 22712, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 53972, 0, 3,
                                                                       49154, 19775, 49484, 5156,
                                                                       5222, 22910, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 54368, 0, 3,
                                                                       49484, 19940, 49814, 5222,
                                                                       5288, 23108, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 54764, 0, 3,
                                                                       49814, 20105, 50144, 5288,
                                                                       5354, 23306, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 55160, 0, 3,
                                                                       50474, 20765, 50804, 5486,
                                                                       5552, 23900, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 55556, 0, 3,
                                                                       50804, 20930, 51134, 5552,
                                                                       5618, 24098, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 55952, 0, 3,
                                                                       51134, 21095, 51464, 5618,
                                                                       5684, 24296, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 56348, 0, 3,
                                                                       51464, 21260, 51794, 5684,
                                                                       5750, 24494, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 56744, 0, 3,
                                                                       51794, 21425, 52124, 5750,
                                                                       5816, 24692, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 57140, 0, 3,
                                                                       52124, 21590, 52454, 5816,
                                                                       5882, 24890, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 57536, 0, 3,
                                                                       52784, 22316, 53180, 6014,
                                                                       6092, 25556, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 58004, 0, 3,
                                                                       53180, 22514, 53576, 6092,
                                                                       6170, 25790, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 58472, 0, 3,
                                                                       53576, 22712, 53972, 6170,
                                                                       6248, 26024, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 58940, 0, 3,
                                                                       53972, 22910, 54368, 6248,
                                                                       6326, 26258, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 59408, 0, 3,
                                                                       54368, 23108, 54764, 6326,
                                                                       6404, 26492, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 59876, 0, 3,
                                                                       55160, 23900, 55556, 6560,
                                                                       6638, 27194, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 60344, 0, 3,
                                                                       55556, 24098, 55952, 6638,
                                                                       6716, 27428, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 60812, 0, 3,
                                                                       55952, 24296, 56348, 6716,
                                                                       6794, 27662, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 61280, 0, 3,
                                                                       56348, 24494, 56744, 6794,
                                                                       6872, 27896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 61748, 0, 3,
                                                                       56744, 24692, 57140, 6872,
                                                                       6950, 28130, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62216, 3, 7106,
                                                                       7109, 28364, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62226, 3, 7109,
                                                                       7112, 28370, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62236, 3, 7112,
                                                                       7115, 28376, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62246, 3, 7115,
                                                                       7118, 28382, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62256, 3, 7118,
                                                                       7121, 28388, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62266, 3, 7121,
                                                                       7124, 28394, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62276, 3, 7124,
                                                                       7127, 28400, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62286, 3, 7127,
                                                                       7130, 28406, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62296, 3, 7130,
                                                                       7133, 28412, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62306, 3, 7133,
                                                                       7136, 28418, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62316, 3, 7136,
                                                                       7139, 28424, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62326, 3, 7139,
                                                                       7142, 28430, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62336, 3, 7142,
                                                                       7145, 28436, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62346, 3, 7145,
                                                                       7148, 28442, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62356, 3, 7148,
                                                                       7151, 28448, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62366, 3, 7151,
                                                                       7154, 28454, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62376, 3, 7160,
                                                                       7163, 28460, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62386, 3, 7163,
                                                                       7166, 28466, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62396, 3, 7166,
                                                                       7169, 28472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62406, 3, 7169,
                                                                       7172, 28478, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62416, 3, 7172,
                                                                       7175, 28484, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62426, 3, 7175,
                                                                       7178, 28490, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62436, 3, 7178,
                                                                       7181, 28496, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62446, 3, 7181,
                                                                       7184, 28502, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62456, 3, 7184,
                                                                       7187, 28508, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62466, 3, 7187,
                                                                       7190, 28514, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62476, 3, 7190,
                                                                       7193, 28520, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62486, 3, 7193,
                                                                       7196, 28526, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62496, 3, 7196,
                                                                       7199, 28532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62506, 3, 7199,
                                                                       7202, 28538, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62516, 3, 7202,
                                                                       7205, 28544, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 62526, 3, 7205,
                                                                       7208, 28550, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62536, 0, 3,
                                                                       62216, 28364, 62226,
                                                                       28556, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62566, 0, 3,
                                                                       62226, 28370, 62236,
                                                                       28574, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62596, 0, 3,
                                                                       62236, 28376, 62246,
                                                                       28592, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62626, 0, 3,
                                                                       62246, 28382, 62256,
                                                                       28610, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62656, 0, 3,
                                                                       62256, 28388, 62266,
                                                                       28628, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62686, 0, 3,
                                                                       62266, 28394, 62276,
                                                                       28646, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62716, 0, 3,
                                                                       62276, 28400, 62286,
                                                                       28664, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62746, 0, 3,
                                                                       62286, 28406, 62296,
                                                                       28682, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62776, 0, 3,
                                                                       62296, 28412, 62306,
                                                                       28700, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62806, 0, 3,
                                                                       62306, 28418, 62316,
                                                                       28718, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62836, 0, 3,
                                                                       62316, 28424, 62326,
                                                                       28736, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62866, 0, 3,
                                                                       62326, 28430, 62336,
                                                                       28754, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62896, 0, 3,
                                                                       62336, 28436, 62346,
                                                                       28772, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62926, 0, 3,
                                                                       62346, 28442, 62356,
                                                                       28790, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62956, 0, 3,
                                                                       62356, 28448, 62366,
                                                                       28808, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 62986, 0, 3,
                                                                       62376, 28460, 62386,
                                                                       28826, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63016, 0, 3,
                                                                       62386, 28466, 62396,
                                                                       28844, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63046, 0, 3,
                                                                       62396, 28472, 62406,
                                                                       28862, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63076, 0, 3,
                                                                       62406, 28478, 62416,
                                                                       28880, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63106, 0, 3,
                                                                       62416, 28484, 62426,
                                                                       28898, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63136, 0, 3,
                                                                       62426, 28490, 62436,
                                                                       28916, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63166, 0, 3,
                                                                       62436, 28496, 62446,
                                                                       28934, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63196, 0, 3,
                                                                       62446, 28502, 62456,
                                                                       28952, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63226, 0, 3,
                                                                       62456, 28508, 62466,
                                                                       28970, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63256, 0, 3,
                                                                       62466, 28514, 62476,
                                                                       28988, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63286, 0, 3,
                                                                       62476, 28520, 62486,
                                                                       29006, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63316, 0, 3,
                                                                       62486, 28526, 62496,
                                                                       29024, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63346, 0, 3,
                                                                       62496, 28532, 62506,
                                                                       29042, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63376, 0, 3,
                                                                       62506, 28538, 62516,
                                                                       29060, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 63406, 0, 3,
                                                                       62516, 28544, 62526,
                                                                       29078, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 63436, 0, 3,
                                                                       62536, 28556, 62566, 7484,
                                                                       7502, 29096, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 63496, 0, 3,
                                                                       62566, 28574, 62596, 7502,
                                                                       7520, 29132, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 63556, 0, 3,
                                                                       62596, 28592, 62626, 7520,
                                                                       7538, 29168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 63616, 0, 3,
                                                                       62626, 28610, 62656, 7538,
                                                                       7556, 29204, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 63676, 0, 3,
                                                                       62656, 28628, 62686, 7556,
                                                                       7574, 29240, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 63736, 0, 3,
                                                                       62686, 28646, 62716, 7574,
                                                                       7592, 29276, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 63796, 0, 3,
                                                                       62716, 28664, 62746, 7592,
                                                                       7610, 29312, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 63856, 0, 3,
                                                                       62746, 28682, 62776, 7610,
                                                                       7628, 29348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 63916, 0, 3,
                                                                       62776, 28700, 62806, 7628,
                                                                       7646, 29384, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 63976, 0, 3,
                                                                       62806, 28718, 62836, 7646,
                                                                       7664, 29420, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64036, 0, 3,
                                                                       62836, 28736, 62866, 7664,
                                                                       7682, 29456, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64096, 0, 3,
                                                                       62866, 28754, 62896, 7682,
                                                                       7700, 29492, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64156, 0, 3,
                                                                       62896, 28772, 62926, 7700,
                                                                       7718, 29528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64216, 0, 3,
                                                                       62926, 28790, 62956, 7718,
                                                                       7736, 29564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64276, 0, 3,
                                                                       62986, 28826, 63016, 7772,
                                                                       7790, 29600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64336, 0, 3,
                                                                       63016, 28844, 63046, 7790,
                                                                       7808, 29636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64396, 0, 3,
                                                                       63046, 28862, 63076, 7808,
                                                                       7826, 29672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64456, 0, 3,
                                                                       63076, 28880, 63106, 7826,
                                                                       7844, 29708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64516, 0, 3,
                                                                       63106, 28898, 63136, 7844,
                                                                       7862, 29744, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64576, 0, 3,
                                                                       63136, 28916, 63166, 7862,
                                                                       7880, 29780, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64636, 0, 3,
                                                                       63166, 28934, 63196, 7880,
                                                                       7898, 29816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64696, 0, 3,
                                                                       63196, 28952, 63226, 7898,
                                                                       7916, 29852, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64756, 0, 3,
                                                                       63226, 28970, 63256, 7916,
                                                                       7934, 29888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64816, 0, 3,
                                                                       63256, 28988, 63286, 7934,
                                                                       7952, 29924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64876, 0, 3,
                                                                       63286, 29006, 63316, 7952,
                                                                       7970, 29960, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64936, 0, 3,
                                                                       63316, 29024, 63346, 7970,
                                                                       7988, 29996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 64996, 0, 3,
                                                                       63346, 29042, 63376, 7988,
                                                                       8006, 30032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 65056, 0, 3,
                                                                       63376, 29060, 63406, 8006,
                                                                       8024, 30068, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 65116, 0, 3,
                                                                       63436, 29096, 63496, 8060,
                                                                       8090, 30104, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 65216, 0, 3,
                                                                       63496, 29132, 63556, 8090,
                                                                       8120, 30164, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 65316, 0, 3,
                                                                       63556, 29168, 63616, 8120,
                                                                       8150, 30224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 65416, 0, 3,
                                                                       63616, 29204, 63676, 8150,
                                                                       8180, 30284, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 65516, 0, 3,
                                                                       63676, 29240, 63736, 8180,
                                                                       8210, 30344, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 65616, 0, 3,
                                                                       63736, 29276, 63796, 8210,
                                                                       8240, 30404, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 65716, 0, 3,
                                                                       63796, 29312, 63856, 8240,
                                                                       8270, 30464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 65816, 0, 3,
                                                                       63856, 29348, 63916, 8270,
                                                                       8300, 30524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 65916, 0, 3,
                                                                       63916, 29384, 63976, 8300,
                                                                       8330, 30584, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 66016, 0, 3,
                                                                       63976, 29420, 64036, 8330,
                                                                       8360, 30644, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 66116, 0, 3,
                                                                       64036, 29456, 64096, 8360,
                                                                       8390, 30704, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 66216, 0, 3,
                                                                       64096, 29492, 64156, 8390,
                                                                       8420, 30764, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 66316, 0, 3,
                                                                       64156, 29528, 64216, 8420,
                                                                       8450, 30824, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 66416, 0, 3,
                                                                       64276, 29600, 64336, 8510,
                                                                       8540, 30884, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 66516, 0, 3,
                                                                       64336, 29636, 64396, 8540,
                                                                       8570, 30944, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 66616, 0, 3,
                                                                       64396, 29672, 64456, 8570,
                                                                       8600, 31004, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 66716, 0, 3,
                                                                       64456, 29708, 64516, 8600,
                                                                       8630, 31064, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 66816, 0, 3,
                                                                       64516, 29744, 64576, 8630,
                                                                       8660, 31124, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 66916, 0, 3,
                                                                       64576, 29780, 64636, 8660,
                                                                       8690, 31184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 67016, 0, 3,
                                                                       64636, 29816, 64696, 8690,
                                                                       8720, 31244, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 67116, 0, 3,
                                                                       64696, 29852, 64756, 8720,
                                                                       8750, 31304, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 67216, 0, 3,
                                                                       64756, 29888, 64816, 8750,
                                                                       8780, 31364, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 67316, 0, 3,
                                                                       64816, 29924, 64876, 8780,
                                                                       8810, 31424, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 67416, 0, 3,
                                                                       64876, 29960, 64936, 8810,
                                                                       8840, 31484, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 67516, 0, 3,
                                                                       64936, 29996, 64996, 8840,
                                                                       8870, 31544, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 67616, 0, 3,
                                                                       64996, 30032, 65056, 8870,
                                                                       8900, 31604, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 67716, 0, 3,
                                                                       65116, 30104, 65216, 8960,
                                                                       9005, 31664, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 67866, 0, 3,
                                                                       65216, 30164, 65316, 9005,
                                                                       9050, 31754, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 68016, 0, 3,
                                                                       65316, 30224, 65416, 9050,
                                                                       9095, 31844, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 68166, 0, 3,
                                                                       65416, 30284, 65516, 9095,
                                                                       9140, 31934, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 68316, 0, 3,
                                                                       65516, 30344, 65616, 9140,
                                                                       9185, 32024, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 68466, 0, 3,
                                                                       65616, 30404, 65716, 9185,
                                                                       9230, 32114, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 68616, 0, 3,
                                                                       65716, 30464, 65816, 9230,
                                                                       9275, 32204, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 68766, 0, 3,
                                                                       65816, 30524, 65916, 9275,
                                                                       9320, 32294, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 68916, 0, 3,
                                                                       65916, 30584, 66016, 9320,
                                                                       9365, 32384, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 69066, 0, 3,
                                                                       66016, 30644, 66116, 9365,
                                                                       9410, 32474, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 69216, 0, 3,
                                                                       66116, 30704, 66216, 9410,
                                                                       9455, 32564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 69366, 0, 3,
                                                                       66216, 30764, 66316, 9455,
                                                                       9500, 32654, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 69516, 0, 3,
                                                                       66416, 30884, 66516, 9590,
                                                                       9635, 32744, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 69666, 0, 3,
                                                                       66516, 30944, 66616, 9635,
                                                                       9680, 32834, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 69816, 0, 3,
                                                                       66616, 31004, 66716, 9680,
                                                                       9725, 32924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 69966, 0, 3,
                                                                       66716, 31064, 66816, 9725,
                                                                       9770, 33014, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 70116, 0, 3,
                                                                       66816, 31124, 66916, 9770,
                                                                       9815, 33104, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 70266, 0, 3,
                                                                       66916, 31184, 67016, 9815,
                                                                       9860, 33194, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 70416, 0, 3,
                                                                       67016, 31244, 67116, 9860,
                                                                       9905, 33284, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 70566, 0, 3,
                                                                       67116, 31304, 67216, 9905,
                                                                       9950, 33374, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 70716, 0, 3,
                                                                       67216, 31364, 67316, 9950,
                                                                       9995, 33464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 70866, 0, 3,
                                                                       67316, 31424, 67416, 9995,
                                                                       10040, 33554, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 71016, 0, 3,
                                                                       67416, 31484, 67516,
                                                                       10040, 10085, 33644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 71166, 0, 3,
                                                                       67516, 31544, 67616,
                                                                       10085, 10130, 33734,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 71316, 0, 3,
                                                                       67716, 31664, 67866,
                                                                       10220, 10283, 33824,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 71526, 0, 3,
                                                                       67866, 31754, 68016,
                                                                       10283, 10346, 33950,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 71736, 0, 3,
                                                                       68016, 31844, 68166,
                                                                       10346, 10409, 34076,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 71946, 0, 3,
                                                                       68166, 31934, 68316,
                                                                       10409, 10472, 34202,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 72156, 0, 3,
                                                                       68316, 32024, 68466,
                                                                       10472, 10535, 34328,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 72366, 0, 3,
                                                                       68466, 32114, 68616,
                                                                       10535, 10598, 34454,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 72576, 0, 3,
                                                                       68616, 32204, 68766,
                                                                       10598, 10661, 34580,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 72786, 0, 3,
                                                                       68766, 32294, 68916,
                                                                       10661, 10724, 34706,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 72996, 0, 3,
                                                                       68916, 32384, 69066,
                                                                       10724, 10787, 34832,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 73206, 0, 3,
                                                                       69066, 32474, 69216,
                                                                       10787, 10850, 34958,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 73416, 0, 3,
                                                                       69216, 32564, 69366,
                                                                       10850, 10913, 35084,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 73626, 0, 3,
                                                                       69516, 32744, 69666,
                                                                       11039, 11102, 35210,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 73836, 0, 3,
                                                                       69666, 32834, 69816,
                                                                       11102, 11165, 35336,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 74046, 0, 3,
                                                                       69816, 32924, 69966,
                                                                       11165, 11228, 35462,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 74256, 0, 3,
                                                                       69966, 33014, 70116,
                                                                       11228, 11291, 35588,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 74466, 0, 3,
                                                                       70116, 33104, 70266,
                                                                       11291, 11354, 35714,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 74676, 0, 3,
                                                                       70266, 33194, 70416,
                                                                       11354, 11417, 35840,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 74886, 0, 3,
                                                                       70416, 33284, 70566,
                                                                       11417, 11480, 35966,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 75096, 0, 3,
                                                                       70566, 33374, 70716,
                                                                       11480, 11543, 36092,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 75306, 0, 3,
                                                                       70716, 33464, 70866,
                                                                       11543, 11606, 36218,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 75516, 0, 3,
                                                                       70866, 33554, 71016,
                                                                       11606, 11669, 36344,
                                                                       ncols, gamma, p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 75726, 0, 3,
                                                                       71016, 33644, 71166,
                                                                       11669, 11732, 36470,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 75936, 0, 3,
                                                                       71316, 33824, 71526,
                                                                       11858, 11942, 36596,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 76216, 0, 3,
                                                                       71526, 33950, 71736,
                                                                       11942, 12026, 36764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 76496, 0, 3,
                                                                       71736, 34076, 71946,
                                                                       12026, 12110, 36932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 76776, 0, 3,
                                                                       71946, 34202, 72156,
                                                                       12110, 12194, 37100,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 77056, 0, 3,
                                                                       72156, 34328, 72366,
                                                                       12194, 12278, 37268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 77336, 0, 3,
                                                                       72366, 34454, 72576,
                                                                       12278, 12362, 37436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 77616, 0, 3,
                                                                       72576, 34580, 72786,
                                                                       12362, 12446, 37604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 77896, 0, 3,
                                                                       72786, 34706, 72996,
                                                                       12446, 12530, 37772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 78176, 0, 3,
                                                                       72996, 34832, 73206,
                                                                       12530, 12614, 37940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 78456, 0, 3,
                                                                       73206, 34958, 73416,
                                                                       12614, 12698, 38108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 78736, 0, 3,
                                                                       73626, 35210, 73836,
                                                                       12866, 12950, 38276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 79016, 0, 3,
                                                                       73836, 35336, 74046,
                                                                       12950, 13034, 38444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 79296, 0, 3,
                                                                       74046, 35462, 74256,
                                                                       13034, 13118, 38612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 79576, 0, 3,
                                                                       74256, 35588, 74466,
                                                                       13118, 13202, 38780,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 79856, 0, 3,
                                                                       74466, 35714, 74676,
                                                                       13202, 13286, 38948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 80136, 0, 3,
                                                                       74676, 35840, 74886,
                                                                       13286, 13370, 39116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 80416, 0, 3,
                                                                       74886, 35966, 75096,
                                                                       13370, 13454, 39284,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 80696, 0, 3,
                                                                       75096, 36092, 75306,
                                                                       13454, 13538, 39452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 80976, 0, 3,
                                                                       75306, 36218, 75516,
                                                                       13538, 13622, 39620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 81256, 0, 3,
                                                                       75516, 36344, 75726,
                                                                       13622, 13706, 39788,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 81536, 0, 3,
                                                                       75936, 36596, 76216,
                                                                       13874, 13982, 39956,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 81896, 0, 3,
                                                                       76216, 36764, 76496,
                                                                       13982, 14090, 40172,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 82256, 0, 3,
                                                                       76496, 36932, 76776,
                                                                       14090, 14198, 40388,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 82616, 0, 3,
                                                                       76776, 37100, 77056,
                                                                       14198, 14306, 40604,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 82976, 0, 3,
                                                                       77056, 37268, 77336,
                                                                       14306, 14414, 40820,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 83336, 0, 3,
                                                                       77336, 37436, 77616,
                                                                       14414, 14522, 41036,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 83696, 0, 3,
                                                                       77616, 37604, 77896,
                                                                       14522, 14630, 41252,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 84056, 0, 3,
                                                                       77896, 37772, 78176,
                                                                       14630, 14738, 41468,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 84416, 0, 3,
                                                                       78176, 37940, 78456,
                                                                       14738, 14846, 41684,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 84776, 0, 3,
                                                                       78736, 38276, 79016,
                                                                       15062, 15170, 41900,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 85136, 0, 3,
                                                                       79016, 38444, 79296,
                                                                       15170, 15278, 42116,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 85496, 0, 3,
                                                                       79296, 38612, 79576,
                                                                       15278, 15386, 42332,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 85856, 0, 3,
                                                                       79576, 38780, 79856,
                                                                       15386, 15494, 42548,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 86216, 0, 3,
                                                                       79856, 38948, 80136,
                                                                       15494, 15602, 42764,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 86576, 0, 3,
                                                                       80136, 39116, 80416,
                                                                       15602, 15710, 42980,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 86936, 0, 3,
                                                                       80416, 39284, 80696,
                                                                       15710, 15818, 43196,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 87296, 0, 3,
                                                                       80696, 39452, 80976,
                                                                       15818, 15926, 43412,
                                                                       ncols, gamma, p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 87656, 0, 3,
                                                                       80976, 39620, 81256,
                                                                       15926, 16034, 43628,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 88016, 0, 3,
                                                                       81536, 39956, 81896,
                                                                       16250, 16385, 43844,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 88466, 0, 3,
                                                                       81896, 40172, 82256,
                                                                       16385, 16520, 44114,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 88916, 0, 3,
                                                                       82256, 40388, 82616,
                                                                       16520, 16655, 44384,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 89366, 0, 3,
                                                                       82616, 40604, 82976,
                                                                       16655, 16790, 44654,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 89816, 0, 3,
                                                                       82976, 40820, 83336,
                                                                       16790, 16925, 44924,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 90266, 0, 3,
                                                                       83336, 41036, 83696,
                                                                       16925, 17060, 45194,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 90716, 0, 3,
                                                                       83696, 41252, 84056,
                                                                       17060, 17195, 45464,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 91166, 0, 3,
                                                                       84056, 41468, 84416,
                                                                       17195, 17330, 45734,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 91616, 0, 3,
                                                                       84776, 41900, 85136,
                                                                       17600, 17735, 46004,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 92066, 0, 3,
                                                                       85136, 42116, 85496,
                                                                       17735, 17870, 46274,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 92516, 0, 3,
                                                                       85496, 42332, 85856,
                                                                       17870, 18005, 46544,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 92966, 0, 3,
                                                                       85856, 42548, 86216,
                                                                       18005, 18140, 46814,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 93416, 0, 3,
                                                                       86216, 42764, 86576,
                                                                       18140, 18275, 47084,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 93866, 0, 3,
                                                                       86576, 42980, 86936,
                                                                       18275, 18410, 47354,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 94316, 0, 3,
                                                                       86936, 43196, 87296,
                                                                       18410, 18545, 47624,
                                                                       ncols, gamma, p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 94766, 0, 3,
                                                                       87296, 43412, 87656,
                                                                       18545, 18680, 47894,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 95216, 0, 3,
                                                                       88016, 43844, 88466,
                                                                       18950, 19115, 48164,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 95766, 0, 3,
                                                                       88466, 44114, 88916,
                                                                       19115, 19280, 48494,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 96316, 0, 3,
                                                                       88916, 44384, 89366,
                                                                       19280, 19445, 48824,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 96866, 0, 3,
                                                                       89366, 44654, 89816,
                                                                       19445, 19610, 49154,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 97416, 0, 3,
                                                                       89816, 44924, 90266,
                                                                       19610, 19775, 49484,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 97966, 0, 3,
                                                                       90266, 45194, 90716,
                                                                       19775, 19940, 49814,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 98516, 0, 3,
                                                                       90716, 45464, 91166,
                                                                       19940, 20105, 50144,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 99066, 0, 3,
                                                                       91616, 46004, 92066,
                                                                       20435, 20600, 50474,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 99616, 0, 3,
                                                                       92066, 46274, 92516,
                                                                       20600, 20765, 50804,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 100166, 0, 3,
                                                                       92516, 46544, 92966,
                                                                       20765, 20930, 51134,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 100716, 0, 3,
                                                                       92966, 46814, 93416,
                                                                       20930, 21095, 51464,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 101266, 0, 3,
                                                                       93416, 47084, 93866,
                                                                       21095, 21260, 51794,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 101816, 0, 3,
                                                                       93866, 47354, 94316,
                                                                       21260, 21425, 52124,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 102366, 0, 3,
                                                                       94316, 47624, 94766,
                                                                       21425, 21590, 52454,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 102916, 0, 3,
                                                                       95216, 48164, 95766,
                                                                       21920, 22118, 52784,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 103576, 0, 3,
                                                                       95766, 48494, 96316,
                                                                       22118, 22316, 53180,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 104236, 0, 3,
                                                                       96316, 48824, 96866,
                                                                       22316, 22514, 53576,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 104896, 0, 3,
                                                                       96866, 49154, 97416,
                                                                       22514, 22712, 53972,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 105556, 0, 3,
                                                                       97416, 49484, 97966,
                                                                       22712, 22910, 54368,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 106216, 0, 3,
                                                                       97966, 49814, 98516,
                                                                       22910, 23108, 54764,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 106876, 0, 3,
                                                                       99066, 50474, 99616,
                                                                       23504, 23702, 55160,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 107536, 0, 3,
                                                                       99616, 50804, 100166,
                                                                       23702, 23900, 55556,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 108196, 0, 3,
                                                                       100166, 51134, 100716,
                                                                       23900, 24098, 55952,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 108856, 0, 3,
                                                                       100716, 51464, 101266,
                                                                       24098, 24296, 56348,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 109516, 0, 3,
                                                                       101266, 51794, 101816,
                                                                       24296, 24494, 56744,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 110176, 0, 3,
                                                                       101816, 52124, 102366,
                                                                       24494, 24692, 57140,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 110836, 0, 3,
                                                                       102916, 52784, 103576,
                                                                       25088, 25322, 57536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 111616, 0, 3,
                                                                       103576, 53180, 104236,
                                                                       25322, 25556, 58004,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 112396, 0, 3,
                                                                       104236, 53576, 104896,
                                                                       25556, 25790, 58472,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 113176, 0, 3,
                                                                       104896, 53972, 105556,
                                                                       25790, 26024, 58940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 113956, 0, 3,
                                                                       105556, 54368, 106216,
                                                                       26024, 26258, 59408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 114736, 0, 3,
                                                                       106876, 55160, 107536,
                                                                       26726, 26960, 59876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 115516, 0, 3,
                                                                       107536, 55556, 108196,
                                                                       26960, 27194, 60344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 116296, 0, 3,
                                                                       108196, 55952, 108856,
                                                                       27194, 27428, 60812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 117076, 0, 3,
                                                                       108856, 56348, 109516,
                                                                       27428, 27662, 61280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 117856, 0, 3,
                                                                       109516, 56744, 110176,
                                                                       27662, 27896, 61748,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118636, 3, 28364,
                                                                       28370, 62236, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118651, 3, 28370,
                                                                       28376, 62246, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118666, 3, 28376,
                                                                       28382, 62256, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118681, 3, 28382,
                                                                       28388, 62266, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118696, 3, 28388,
                                                                       28394, 62276, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118711, 3, 28394,
                                                                       28400, 62286, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118726, 3, 28400,
                                                                       28406, 62296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118741, 3, 28406,
                                                                       28412, 62306, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118756, 3, 28412,
                                                                       28418, 62316, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118771, 3, 28418,
                                                                       28424, 62326, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118786, 3, 28424,
                                                                       28430, 62336, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118801, 3, 28430,
                                                                       28436, 62346, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118816, 3, 28436,
                                                                       28442, 62356, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118831, 3, 28442,
                                                                       28448, 62366, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118846, 3, 28460,
                                                                       28466, 62396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118861, 3, 28466,
                                                                       28472, 62406, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118876, 3, 28472,
                                                                       28478, 62416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118891, 3, 28478,
                                                                       28484, 62426, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118906, 3, 28484,
                                                                       28490, 62436, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118921, 3, 28490,
                                                                       28496, 62446, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118936, 3, 28496,
                                                                       28502, 62456, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118951, 3, 28502,
                                                                       28508, 62466, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118966, 3, 28508,
                                                                       28514, 62476, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118981, 3, 28514,
                                                                       28520, 62486, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 118996, 3, 28520,
                                                                       28526, 62496, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 119011, 3, 28526,
                                                                       28532, 62506, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 119026, 3, 28532,
                                                                       28538, 62516, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 119041, 3, 28538,
                                                                       28544, 62526, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119056, 0, 3,
                                                                       118636, 62236, 118651,
                                                                       28556, 28574, 62596,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119101, 0, 3,
                                                                       118651, 62246, 118666,
                                                                       28574, 28592, 62626,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119146, 0, 3,
                                                                       118666, 62256, 118681,
                                                                       28592, 28610, 62656,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119191, 0, 3,
                                                                       118681, 62266, 118696,
                                                                       28610, 28628, 62686,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119236, 0, 3,
                                                                       118696, 62276, 118711,
                                                                       28628, 28646, 62716,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119281, 0, 3,
                                                                       118711, 62286, 118726,
                                                                       28646, 28664, 62746,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119326, 0, 3,
                                                                       118726, 62296, 118741,
                                                                       28664, 28682, 62776,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119371, 0, 3,
                                                                       118741, 62306, 118756,
                                                                       28682, 28700, 62806,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119416, 0, 3,
                                                                       118756, 62316, 118771,
                                                                       28700, 28718, 62836,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119461, 0, 3,
                                                                       118771, 62326, 118786,
                                                                       28718, 28736, 62866,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119506, 0, 3,
                                                                       118786, 62336, 118801,
                                                                       28736, 28754, 62896,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119551, 0, 3,
                                                                       118801, 62346, 118816,
                                                                       28754, 28772, 62926,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119596, 0, 3,
                                                                       118816, 62356, 118831,
                                                                       28772, 28790, 62956,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119641, 0, 3,
                                                                       118846, 62396, 118861,
                                                                       28826, 28844, 63046,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119686, 0, 3,
                                                                       118861, 62406, 118876,
                                                                       28844, 28862, 63076,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119731, 0, 3,
                                                                       118876, 62416, 118891,
                                                                       28862, 28880, 63106,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119776, 0, 3,
                                                                       118891, 62426, 118906,
                                                                       28880, 28898, 63136,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119821, 0, 3,
                                                                       118906, 62436, 118921,
                                                                       28898, 28916, 63166,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119866, 0, 3,
                                                                       118921, 62446, 118936,
                                                                       28916, 28934, 63196,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119911, 0, 3,
                                                                       118936, 62456, 118951,
                                                                       28934, 28952, 63226,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 119956, 0, 3,
                                                                       118951, 62466, 118966,
                                                                       28952, 28970, 63256,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 120001, 0, 3,
                                                                       118966, 62476, 118981,
                                                                       28970, 28988, 63286,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 120046, 0, 3,
                                                                       118981, 62486, 118996,
                                                                       28988, 29006, 63316,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 120091, 0, 3,
                                                                       118996, 62496, 119011,
                                                                       29006, 29024, 63346,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 120136, 0, 3,
                                                                       119011, 62506, 119026,
                                                                       29024, 29042, 63376,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 120181, 0, 3,
                                                                       119026, 62516, 119041,
                                                                       29042, 29060, 63406,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 120226, 0, 3,
                                                                       119056, 62596, 119101,
                                                                       29096, 29132, 63556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 120316, 0, 3,
                                                                       119101, 62626, 119146,
                                                                       29132, 29168, 63616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 120406, 0, 3,
                                                                       119146, 62656, 119191,
                                                                       29168, 29204, 63676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 120496, 0, 3,
                                                                       119191, 62686, 119236,
                                                                       29204, 29240, 63736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 120586, 0, 3,
                                                                       119236, 62716, 119281,
                                                                       29240, 29276, 63796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 120676, 0, 3,
                                                                       119281, 62746, 119326,
                                                                       29276, 29312, 63856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 120766, 0, 3,
                                                                       119326, 62776, 119371,
                                                                       29312, 29348, 63916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 120856, 0, 3,
                                                                       119371, 62806, 119416,
                                                                       29348, 29384, 63976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 120946, 0, 3,
                                                                       119416, 62836, 119461,
                                                                       29384, 29420, 64036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121036, 0, 3,
                                                                       119461, 62866, 119506,
                                                                       29420, 29456, 64096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121126, 0, 3,
                                                                       119506, 62896, 119551,
                                                                       29456, 29492, 64156,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121216, 0, 3,
                                                                       119551, 62926, 119596,
                                                                       29492, 29528, 64216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121306, 0, 3,
                                                                       119641, 63046, 119686,
                                                                       29600, 29636, 64396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121396, 0, 3,
                                                                       119686, 63076, 119731,
                                                                       29636, 29672, 64456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121486, 0, 3,
                                                                       119731, 63106, 119776,
                                                                       29672, 29708, 64516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121576, 0, 3,
                                                                       119776, 63136, 119821,
                                                                       29708, 29744, 64576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121666, 0, 3,
                                                                       119821, 63166, 119866,
                                                                       29744, 29780, 64636,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121756, 0, 3,
                                                                       119866, 63196, 119911,
                                                                       29780, 29816, 64696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121846, 0, 3,
                                                                       119911, 63226, 119956,
                                                                       29816, 29852, 64756,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 121936, 0, 3,
                                                                       119956, 63256, 120001,
                                                                       29852, 29888, 64816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 122026, 0, 3,
                                                                       120001, 63286, 120046,
                                                                       29888, 29924, 64876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 122116, 0, 3,
                                                                       120046, 63316, 120091,
                                                                       29924, 29960, 64936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 122206, 0, 3,
                                                                       120091, 63346, 120136,
                                                                       29960, 29996, 64996,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 122296, 0, 3,
                                                                       120136, 63376, 120181,
                                                                       29996, 30032, 65056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 122386, 0, 3,
                                                                       120226, 63556, 120316,
                                                                       30104, 30164, 65316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 122536, 0, 3,
                                                                       120316, 63616, 120406,
                                                                       30164, 30224, 65416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 122686, 0, 3,
                                                                       120406, 63676, 120496,
                                                                       30224, 30284, 65516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 122836, 0, 3,
                                                                       120496, 63736, 120586,
                                                                       30284, 30344, 65616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 122986, 0, 3,
                                                                       120586, 63796, 120676,
                                                                       30344, 30404, 65716,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 123136, 0, 3,
                                                                       120676, 63856, 120766,
                                                                       30404, 30464, 65816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 123286, 0, 3,
                                                                       120766, 63916, 120856,
                                                                       30464, 30524, 65916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 123436, 0, 3,
                                                                       120856, 63976, 120946,
                                                                       30524, 30584, 66016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 123586, 0, 3,
                                                                       120946, 64036, 121036,
                                                                       30584, 30644, 66116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 123736, 0, 3,
                                                                       121036, 64096, 121126,
                                                                       30644, 30704, 66216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 123886, 0, 3,
                                                                       121126, 64156, 121216,
                                                                       30704, 30764, 66316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 124036, 0, 3,
                                                                       121306, 64396, 121396,
                                                                       30884, 30944, 66616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 124186, 0, 3,
                                                                       121396, 64456, 121486,
                                                                       30944, 31004, 66716,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 124336, 0, 3,
                                                                       121486, 64516, 121576,
                                                                       31004, 31064, 66816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 124486, 0, 3,
                                                                       121576, 64576, 121666,
                                                                       31064, 31124, 66916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 124636, 0, 3,
                                                                       121666, 64636, 121756,
                                                                       31124, 31184, 67016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 124786, 0, 3,
                                                                       121756, 64696, 121846,
                                                                       31184, 31244, 67116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 124936, 0, 3,
                                                                       121846, 64756, 121936,
                                                                       31244, 31304, 67216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 125086, 0, 3,
                                                                       121936, 64816, 122026,
                                                                       31304, 31364, 67316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 125236, 0, 3,
                                                                       122026, 64876, 122116,
                                                                       31364, 31424, 67416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 125386, 0, 3,
                                                                       122116, 64936, 122206,
                                                                       31424, 31484, 67516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 125536, 0, 3,
                                                                       122206, 64996, 122296,
                                                                       31484, 31544, 67616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 125686, 0, 3,
                                                                       122386, 65316, 122536,
                                                                       31664, 31754, 68016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 125911, 0, 3,
                                                                       122536, 65416, 122686,
                                                                       31754, 31844, 68166,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 126136, 0, 3,
                                                                       122686, 65516, 122836,
                                                                       31844, 31934, 68316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 126361, 0, 3,
                                                                       122836, 65616, 122986,
                                                                       31934, 32024, 68466,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 126586, 0, 3,
                                                                       122986, 65716, 123136,
                                                                       32024, 32114, 68616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 126811, 0, 3,
                                                                       123136, 65816, 123286,
                                                                       32114, 32204, 68766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 127036, 0, 3,
                                                                       123286, 65916, 123436,
                                                                       32204, 32294, 68916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 127261, 0, 3,
                                                                       123436, 66016, 123586,
                                                                       32294, 32384, 69066,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 127486, 0, 3,
                                                                       123586, 66116, 123736,
                                                                       32384, 32474, 69216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 127711, 0, 3,
                                                                       123736, 66216, 123886,
                                                                       32474, 32564, 69366,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 127936, 0, 3,
                                                                       124036, 66616, 124186,
                                                                       32744, 32834, 69816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 128161, 0, 3,
                                                                       124186, 66716, 124336,
                                                                       32834, 32924, 69966,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 128386, 0, 3,
                                                                       124336, 66816, 124486,
                                                                       32924, 33014, 70116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 128611, 0, 3,
                                                                       124486, 66916, 124636,
                                                                       33014, 33104, 70266,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 128836, 0, 3,
                                                                       124636, 67016, 124786,
                                                                       33104, 33194, 70416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 129061, 0, 3,
                                                                       124786, 67116, 124936,
                                                                       33194, 33284, 70566,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 129286, 0, 3,
                                                                       124936, 67216, 125086,
                                                                       33284, 33374, 70716,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 129511, 0, 3,
                                                                       125086, 67316, 125236,
                                                                       33374, 33464, 70866,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 129736, 0, 3,
                                                                       125236, 67416, 125386,
                                                                       33464, 33554, 71016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 129961, 0, 3,
                                                                       125386, 67516, 125536,
                                                                       33554, 33644, 71166,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 130186, 0, 3,
                                                                       125686, 68016, 125911,
                                                                       33824, 33950, 71736,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 130501, 0, 3,
                                                                       125911, 68166, 126136,
                                                                       33950, 34076, 71946,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 130816, 0, 3,
                                                                       126136, 68316, 126361,
                                                                       34076, 34202, 72156,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 131131, 0, 3,
                                                                       126361, 68466, 126586,
                                                                       34202, 34328, 72366,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 131446, 0, 3,
                                                                       126586, 68616, 126811,
                                                                       34328, 34454, 72576,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 131761, 0, 3,
                                                                       126811, 68766, 127036,
                                                                       34454, 34580, 72786,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 132076, 0, 3,
                                                                       127036, 68916, 127261,
                                                                       34580, 34706, 72996,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 132391, 0, 3,
                                                                       127261, 69066, 127486,
                                                                       34706, 34832, 73206,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 132706, 0, 3,
                                                                       127486, 69216, 127711,
                                                                       34832, 34958, 73416,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 133021, 0, 3,
                                                                       127936, 69816, 128161,
                                                                       35210, 35336, 74046,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 133336, 0, 3,
                                                                       128161, 69966, 128386,
                                                                       35336, 35462, 74256,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 133651, 0, 3,
                                                                       128386, 70116, 128611,
                                                                       35462, 35588, 74466,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 133966, 0, 3,
                                                                       128611, 70266, 128836,
                                                                       35588, 35714, 74676,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 134281, 0, 3,
                                                                       128836, 70416, 129061,
                                                                       35714, 35840, 74886,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 134596, 0, 3,
                                                                       129061, 70566, 129286,
                                                                       35840, 35966, 75096,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 134911, 0, 3,
                                                                       129286, 70716, 129511,
                                                                       35966, 36092, 75306,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 135226, 0, 3,
                                                                       129511, 70866, 129736,
                                                                       36092, 36218, 75516,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 135541, 0, 3,
                                                                       129736, 71016, 129961,
                                                                       36218, 36344, 75726,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 135856, 0, 3,
                                                                       130186, 71736, 130501,
                                                                       36596, 36764, 76496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 136276, 0, 3,
                                                                       130501, 71946, 130816,
                                                                       36764, 36932, 76776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 136696, 0, 3,
                                                                       130816, 72156, 131131,
                                                                       36932, 37100, 77056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 137116, 0, 3,
                                                                       131131, 72366, 131446,
                                                                       37100, 37268, 77336,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 137536, 0, 3,
                                                                       131446, 72576, 131761,
                                                                       37268, 37436, 77616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 137956, 0, 3,
                                                                       131761, 72786, 132076,
                                                                       37436, 37604, 77896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 138376, 0, 3,
                                                                       132076, 72996, 132391,
                                                                       37604, 37772, 78176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 138796, 0, 3,
                                                                       132391, 73206, 132706,
                                                                       37772, 37940, 78456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 139216, 0, 3,
                                                                       133021, 74046, 133336,
                                                                       38276, 38444, 79296,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 139636, 0, 3,
                                                                       133336, 74256, 133651,
                                                                       38444, 38612, 79576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 140056, 0, 3,
                                                                       133651, 74466, 133966,
                                                                       38612, 38780, 79856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 140476, 0, 3,
                                                                       133966, 74676, 134281,
                                                                       38780, 38948, 80136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 140896, 0, 3,
                                                                       134281, 74886, 134596,
                                                                       38948, 39116, 80416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 141316, 0, 3,
                                                                       134596, 75096, 134911,
                                                                       39116, 39284, 80696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 141736, 0, 3,
                                                                       134911, 75306, 135226,
                                                                       39284, 39452, 80976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 142156, 0, 3,
                                                                       135226, 75516, 135541,
                                                                       39452, 39620, 81256,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 142576, 0, 3,
                                                                       135856, 76496, 136276,
                                                                       39956, 40172, 82256,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 143116, 0, 3,
                                                                       136276, 76776, 136696,
                                                                       40172, 40388, 82616,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 143656, 0, 3,
                                                                       136696, 77056, 137116,
                                                                       40388, 40604, 82976,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 144196, 0, 3,
                                                                       137116, 77336, 137536,
                                                                       40604, 40820, 83336,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 144736, 0, 3,
                                                                       137536, 77616, 137956,
                                                                       40820, 41036, 83696,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 145276, 0, 3,
                                                                       137956, 77896, 138376,
                                                                       41036, 41252, 84056,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 145816, 0, 3,
                                                                       138376, 78176, 138796,
                                                                       41252, 41468, 84416,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 146356, 0, 3,
                                                                       139216, 79296, 139636,
                                                                       41900, 42116, 85496,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 146896, 0, 3,
                                                                       139636, 79576, 140056,
                                                                       42116, 42332, 85856,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 147436, 0, 3,
                                                                       140056, 79856, 140476,
                                                                       42332, 42548, 86216,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 147976, 0, 3,
                                                                       140476, 80136, 140896,
                                                                       42548, 42764, 86576,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 148516, 0, 3,
                                                                       140896, 80416, 141316,
                                                                       42764, 42980, 86936,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 149056, 0, 3,
                                                                       141316, 80696, 141736,
                                                                       42980, 43196, 87296,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 149596, 0, 3,
                                                                       141736, 80976, 142156,
                                                                       43196, 43412, 87656,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 150136, 0, 3,
                                                                       142576, 82256, 143116,
                                                                       43844, 44114, 88916,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 150811, 0, 3,
                                                                       143116, 82616, 143656,
                                                                       44114, 44384, 89366,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 151486, 0, 3,
                                                                       143656, 82976, 144196,
                                                                       44384, 44654, 89816,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 152161, 0, 3,
                                                                       144196, 83336, 144736,
                                                                       44654, 44924, 90266,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 152836, 0, 3,
                                                                       144736, 83696, 145276,
                                                                       44924, 45194, 90716,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 153511, 0, 3,
                                                                       145276, 84056, 145816,
                                                                       45194, 45464, 91166,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 154186, 0, 3,
                                                                       146356, 85496, 146896,
                                                                       46004, 46274, 92516,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 154861, 0, 3,
                                                                       146896, 85856, 147436,
                                                                       46274, 46544, 92966,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 155536, 0, 3,
                                                                       147436, 86216, 147976,
                                                                       46544, 46814, 93416,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 156211, 0, 3,
                                                                       147976, 86576, 148516,
                                                                       46814, 47084, 93866,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 156886, 0, 3,
                                                                       148516, 86936, 149056,
                                                                       47084, 47354, 94316,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 157561, 0, 3,
                                                                       149056, 87296, 149596,
                                                                       47354, 47624, 94766,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 158236, 0, 3,
                                                                       150136, 88916, 150811,
                                                                       48164, 48494, 96316,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 159061, 0, 3,
                                                                       150811, 89366, 151486,
                                                                       48494, 48824, 96866,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 159886, 0, 3,
                                                                       151486, 89816, 152161,
                                                                       48824, 49154, 97416,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 160711, 0, 3,
                                                                       152161, 90266, 152836,
                                                                       49154, 49484, 97966,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 161536, 0, 3,
                                                                       152836, 90716, 153511,
                                                                       49484, 49814, 98516,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 162361, 0, 3,
                                                                       154186, 92516, 154861,
                                                                       50474, 50804, 100166,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 163186, 0, 3,
                                                                       154861, 92966, 155536,
                                                                       50804, 51134, 100716,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 164011, 0, 3,
                                                                       155536, 93416, 156211,
                                                                       51134, 51464, 101266,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 164836, 0, 3,
                                                                       156211, 93866, 156886,
                                                                       51464, 51794, 101816,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 165661, 0, 3,
                                                                       156886, 94316, 157561,
                                                                       51794, 52124, 102366,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 166486, 0, 3,
                                                                       158236, 96316, 159061,
                                                                       52784, 53180, 104236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 167476, 0, 3,
                                                                       159061, 96866, 159886,
                                                                       53180, 53576, 104896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 168466, 0, 3,
                                                                       159886, 97416, 160711,
                                                                       53576, 53972, 105556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 169456, 0, 3,
                                                                       160711, 97966, 161536,
                                                                       53972, 54368, 106216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 170446, 0, 3,
                                                                       162361, 100166, 163186,
                                                                       55160, 55556, 108196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 171436, 0, 3,
                                                                       163186, 100716, 164011,
                                                                       55556, 55952, 108856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 172426, 0, 3,
                                                                       164011, 101266, 164836,
                                                                       55952, 56348, 109516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 173416, 0, 3,
                                                                       164836, 101816, 165661,
                                                                       56348, 56744, 110176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 174406, 0, 3,
                                                                       166486, 104236, 167476,
                                                                       57536, 58004, 112396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 175576, 0, 3,
                                                                       167476, 104896, 168466,
                                                                       58004, 58472, 113176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 176746, 0, 3,
                                                                       168466, 105556, 169456,
                                                                       58472, 58940, 113956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 177916, 0, 3,
                                                                       170446, 108196, 171436,
                                                                       59876, 60344, 116296,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 179086, 0, 3,
                                                                       171436, 108856, 172426,
                                                                       60344, 60812, 117076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 180256, 0, 3,
                                                                       172426, 109516, 173416,
                                                                       60812, 61280, 117856,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181426, 3, 62216,
                                                                       62226, 118636, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181447, 3, 62226,
                                                                       62236, 118651, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181468, 3, 62236,
                                                                       62246, 118666, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181489, 3, 62246,
                                                                       62256, 118681, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181510, 3, 62256,
                                                                       62266, 118696, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181531, 3, 62266,
                                                                       62276, 118711, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181552, 3, 62276,
                                                                       62286, 118726, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181573, 3, 62286,
                                                                       62296, 118741, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181594, 3, 62296,
                                                                       62306, 118756, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181615, 3, 62306,
                                                                       62316, 118771, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181636, 3, 62316,
                                                                       62326, 118786, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181657, 3, 62326,
                                                                       62336, 118801, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181678, 3, 62336,
                                                                       62346, 118816, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181699, 3, 62346,
                                                                       62356, 118831, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181720, 3, 62376,
                                                                       62386, 118846, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181741, 3, 62386,
                                                                       62396, 118861, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181762, 3, 62396,
                                                                       62406, 118876, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181783, 3, 62406,
                                                                       62416, 118891, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181804, 3, 62416,
                                                                       62426, 118906, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181825, 3, 62426,
                                                                       62436, 118921, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181846, 3, 62436,
                                                                       62446, 118936, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181867, 3, 62446,
                                                                       62456, 118951, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181888, 3, 62456,
                                                                       62466, 118966, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181909, 3, 62466,
                                                                       62476, 118981, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181930, 3, 62476,
                                                                       62486, 118996, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181951, 3, 62486,
                                                                       62496, 119011, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181972, 3, 62496,
                                                                       62506, 119026, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 181993, 3, 62506,
                                                                       62516, 119041, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182014, 0, 3,
                                                                       181426, 118636, 181447,
                                                                       62536, 62566, 119056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182077, 0, 3,
                                                                       181447, 118651, 181468,
                                                                       62566, 62596, 119101,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182140, 0, 3,
                                                                       181468, 118666, 181489,
                                                                       62596, 62626, 119146,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182203, 0, 3,
                                                                       181489, 118681, 181510,
                                                                       62626, 62656, 119191,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182266, 0, 3,
                                                                       181510, 118696, 181531,
                                                                       62656, 62686, 119236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182329, 0, 3,
                                                                       181531, 118711, 181552,
                                                                       62686, 62716, 119281,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182392, 0, 3,
                                                                       181552, 118726, 181573,
                                                                       62716, 62746, 119326,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182455, 0, 3,
                                                                       181573, 118741, 181594,
                                                                       62746, 62776, 119371,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182518, 0, 3,
                                                                       181594, 118756, 181615,
                                                                       62776, 62806, 119416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182581, 0, 3,
                                                                       181615, 118771, 181636,
                                                                       62806, 62836, 119461,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182644, 0, 3,
                                                                       181636, 118786, 181657,
                                                                       62836, 62866, 119506,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182707, 0, 3,
                                                                       181657, 118801, 181678,
                                                                       62866, 62896, 119551,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182770, 0, 3,
                                                                       181678, 118816, 181699,
                                                                       62896, 62926, 119596,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182833, 0, 3,
                                                                       181720, 118846, 181741,
                                                                       62986, 63016, 119641,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182896, 0, 3,
                                                                       181741, 118861, 181762,
                                                                       63016, 63046, 119686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 182959, 0, 3,
                                                                       181762, 118876, 181783,
                                                                       63046, 63076, 119731,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 183022, 0, 3,
                                                                       181783, 118891, 181804,
                                                                       63076, 63106, 119776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 183085, 0, 3,
                                                                       181804, 118906, 181825,
                                                                       63106, 63136, 119821,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 183148, 0, 3,
                                                                       181825, 118921, 181846,
                                                                       63136, 63166, 119866,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 183211, 0, 3,
                                                                       181846, 118936, 181867,
                                                                       63166, 63196, 119911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 183274, 0, 3,
                                                                       181867, 118951, 181888,
                                                                       63196, 63226, 119956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 183337, 0, 3,
                                                                       181888, 118966, 181909,
                                                                       63226, 63256, 120001,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 183400, 0, 3,
                                                                       181909, 118981, 181930,
                                                                       63256, 63286, 120046,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 183463, 0, 3,
                                                                       181930, 118996, 181951,
                                                                       63286, 63316, 120091,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 183526, 0, 3,
                                                                       181951, 119011, 181972,
                                                                       63316, 63346, 120136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 183589, 0, 3,
                                                                       181972, 119026, 181993,
                                                                       63346, 63376, 120181,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 183652, 0, 3,
                                                                       182014, 119056, 182077,
                                                                       63436, 63496, 120226,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 183778, 0, 3,
                                                                       182077, 119101, 182140,
                                                                       63496, 63556, 120316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 183904, 0, 3,
                                                                       182140, 119146, 182203,
                                                                       63556, 63616, 120406,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 184030, 0, 3,
                                                                       182203, 119191, 182266,
                                                                       63616, 63676, 120496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 184156, 0, 3,
                                                                       182266, 119236, 182329,
                                                                       63676, 63736, 120586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 184282, 0, 3,
                                                                       182329, 119281, 182392,
                                                                       63736, 63796, 120676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 184408, 0, 3,
                                                                       182392, 119326, 182455,
                                                                       63796, 63856, 120766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 184534, 0, 3,
                                                                       182455, 119371, 182518,
                                                                       63856, 63916, 120856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 184660, 0, 3,
                                                                       182518, 119416, 182581,
                                                                       63916, 63976, 120946,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 184786, 0, 3,
                                                                       182581, 119461, 182644,
                                                                       63976, 64036, 121036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 184912, 0, 3,
                                                                       182644, 119506, 182707,
                                                                       64036, 64096, 121126,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 185038, 0, 3,
                                                                       182707, 119551, 182770,
                                                                       64096, 64156, 121216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 185164, 0, 3,
                                                                       182833, 119641, 182896,
                                                                       64276, 64336, 121306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 185290, 0, 3,
                                                                       182896, 119686, 182959,
                                                                       64336, 64396, 121396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 185416, 0, 3,
                                                                       182959, 119731, 183022,
                                                                       64396, 64456, 121486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 185542, 0, 3,
                                                                       183022, 119776, 183085,
                                                                       64456, 64516, 121576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 185668, 0, 3,
                                                                       183085, 119821, 183148,
                                                                       64516, 64576, 121666,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 185794, 0, 3,
                                                                       183148, 119866, 183211,
                                                                       64576, 64636, 121756,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 185920, 0, 3,
                                                                       183211, 119911, 183274,
                                                                       64636, 64696, 121846,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 186046, 0, 3,
                                                                       183274, 119956, 183337,
                                                                       64696, 64756, 121936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 186172, 0, 3,
                                                                       183337, 120001, 183400,
                                                                       64756, 64816, 122026,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 186298, 0, 3,
                                                                       183400, 120046, 183463,
                                                                       64816, 64876, 122116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 186424, 0, 3,
                                                                       183463, 120091, 183526,
                                                                       64876, 64936, 122206,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 186550, 0, 3,
                                                                       183526, 120136, 183589,
                                                                       64936, 64996, 122296,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 186676, 0, 3,
                                                                       183652, 120226, 183778,
                                                                       65116, 65216, 122386,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 186886, 0, 3,
                                                                       183778, 120316, 183904,
                                                                       65216, 65316, 122536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 187096, 0, 3,
                                                                       183904, 120406, 184030,
                                                                       65316, 65416, 122686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 187306, 0, 3,
                                                                       184030, 120496, 184156,
                                                                       65416, 65516, 122836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 187516, 0, 3,
                                                                       184156, 120586, 184282,
                                                                       65516, 65616, 122986,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 187726, 0, 3,
                                                                       184282, 120676, 184408,
                                                                       65616, 65716, 123136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 187936, 0, 3,
                                                                       184408, 120766, 184534,
                                                                       65716, 65816, 123286,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 188146, 0, 3,
                                                                       184534, 120856, 184660,
                                                                       65816, 65916, 123436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 188356, 0, 3,
                                                                       184660, 120946, 184786,
                                                                       65916, 66016, 123586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 188566, 0, 3,
                                                                       184786, 121036, 184912,
                                                                       66016, 66116, 123736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 188776, 0, 3,
                                                                       184912, 121126, 185038,
                                                                       66116, 66216, 123886,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 188986, 0, 3,
                                                                       185164, 121306, 185290,
                                                                       66416, 66516, 124036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 189196, 0, 3,
                                                                       185290, 121396, 185416,
                                                                       66516, 66616, 124186,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 189406, 0, 3,
                                                                       185416, 121486, 185542,
                                                                       66616, 66716, 124336,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 189616, 0, 3,
                                                                       185542, 121576, 185668,
                                                                       66716, 66816, 124486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 189826, 0, 3,
                                                                       185668, 121666, 185794,
                                                                       66816, 66916, 124636,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 190036, 0, 3,
                                                                       185794, 121756, 185920,
                                                                       66916, 67016, 124786,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 190246, 0, 3,
                                                                       185920, 121846, 186046,
                                                                       67016, 67116, 124936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 190456, 0, 3,
                                                                       186046, 121936, 186172,
                                                                       67116, 67216, 125086,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 190666, 0, 3,
                                                                       186172, 122026, 186298,
                                                                       67216, 67316, 125236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 190876, 0, 3,
                                                                       186298, 122116, 186424,
                                                                       67316, 67416, 125386,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 191086, 0, 3,
                                                                       186424, 122206, 186550,
                                                                       67416, 67516, 125536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 191296, 0, 3,
                                                                       186676, 122386, 186886,
                                                                       67716, 67866, 125686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 191611, 0, 3,
                                                                       186886, 122536, 187096,
                                                                       67866, 68016, 125911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 191926, 0, 3,
                                                                       187096, 122686, 187306,
                                                                       68016, 68166, 126136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 192241, 0, 3,
                                                                       187306, 122836, 187516,
                                                                       68166, 68316, 126361,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 192556, 0, 3,
                                                                       187516, 122986, 187726,
                                                                       68316, 68466, 126586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 192871, 0, 3,
                                                                       187726, 123136, 187936,
                                                                       68466, 68616, 126811,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 193186, 0, 3,
                                                                       187936, 123286, 188146,
                                                                       68616, 68766, 127036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 193501, 0, 3,
                                                                       188146, 123436, 188356,
                                                                       68766, 68916, 127261,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 193816, 0, 3,
                                                                       188356, 123586, 188566,
                                                                       68916, 69066, 127486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 194131, 0, 3,
                                                                       188566, 123736, 188776,
                                                                       69066, 69216, 127711,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 194446, 0, 3,
                                                                       188986, 124036, 189196,
                                                                       69516, 69666, 127936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 194761, 0, 3,
                                                                       189196, 124186, 189406,
                                                                       69666, 69816, 128161,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 195076, 0, 3,
                                                                       189406, 124336, 189616,
                                                                       69816, 69966, 128386,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 195391, 0, 3,
                                                                       189616, 124486, 189826,
                                                                       69966, 70116, 128611,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 195706, 0, 3,
                                                                       189826, 124636, 190036,
                                                                       70116, 70266, 128836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 196021, 0, 3,
                                                                       190036, 124786, 190246,
                                                                       70266, 70416, 129061,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 196336, 0, 3,
                                                                       190246, 124936, 190456,
                                                                       70416, 70566, 129286,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 196651, 0, 3,
                                                                       190456, 125086, 190666,
                                                                       70566, 70716, 129511,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 196966, 0, 3,
                                                                       190666, 125236, 190876,
                                                                       70716, 70866, 129736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 197281, 0, 3,
                                                                       190876, 125386, 191086,
                                                                       70866, 71016, 129961,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 197596, 0, 3,
                                                                       191296, 125686, 191611,
                                                                       71316, 71526, 130186,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 198037, 0, 3,
                                                                       191611, 125911, 191926,
                                                                       71526, 71736, 130501,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 198478, 0, 3,
                                                                       191926, 126136, 192241,
                                                                       71736, 71946, 130816,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 198919, 0, 3,
                                                                       192241, 126361, 192556,
                                                                       71946, 72156, 131131,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 199360, 0, 3,
                                                                       192556, 126586, 192871,
                                                                       72156, 72366, 131446,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 199801, 0, 3,
                                                                       192871, 126811, 193186,
                                                                       72366, 72576, 131761,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 200242, 0, 3,
                                                                       193186, 127036, 193501,
                                                                       72576, 72786, 132076,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 200683, 0, 3,
                                                                       193501, 127261, 193816,
                                                                       72786, 72996, 132391,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 201124, 0, 3,
                                                                       193816, 127486, 194131,
                                                                       72996, 73206, 132706,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 201565, 0, 3,
                                                                       194446, 127936, 194761,
                                                                       73626, 73836, 133021,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 202006, 0, 3,
                                                                       194761, 128161, 195076,
                                                                       73836, 74046, 133336,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 202447, 0, 3,
                                                                       195076, 128386, 195391,
                                                                       74046, 74256, 133651,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 202888, 0, 3,
                                                                       195391, 128611, 195706,
                                                                       74256, 74466, 133966,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 203329, 0, 3,
                                                                       195706, 128836, 196021,
                                                                       74466, 74676, 134281,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 203770, 0, 3,
                                                                       196021, 129061, 196336,
                                                                       74676, 74886, 134596,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 204211, 0, 3,
                                                                       196336, 129286, 196651,
                                                                       74886, 75096, 134911,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 204652, 0, 3,
                                                                       196651, 129511, 196966,
                                                                       75096, 75306, 135226,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 205093, 0, 3,
                                                                       196966, 129736, 197281,
                                                                       75306, 75516, 135541,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 205534, 0, 3,
                                                                       197596, 130186, 198037,
                                                                       75936, 76216, 135856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 206122, 0, 3,
                                                                       198037, 130501, 198478,
                                                                       76216, 76496, 136276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 206710, 0, 3,
                                                                       198478, 130816, 198919,
                                                                       76496, 76776, 136696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 207298, 0, 3,
                                                                       198919, 131131, 199360,
                                                                       76776, 77056, 137116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 207886, 0, 3,
                                                                       199360, 131446, 199801,
                                                                       77056, 77336, 137536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 208474, 0, 3,
                                                                       199801, 131761, 200242,
                                                                       77336, 77616, 137956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 209062, 0, 3,
                                                                       200242, 132076, 200683,
                                                                       77616, 77896, 138376,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 209650, 0, 3,
                                                                       200683, 132391, 201124,
                                                                       77896, 78176, 138796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 210238, 0, 3,
                                                                       201565, 133021, 202006,
                                                                       78736, 79016, 139216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 210826, 0, 3,
                                                                       202006, 133336, 202447,
                                                                       79016, 79296, 139636,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 211414, 0, 3,
                                                                       202447, 133651, 202888,
                                                                       79296, 79576, 140056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 212002, 0, 3,
                                                                       202888, 133966, 203329,
                                                                       79576, 79856, 140476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 212590, 0, 3,
                                                                       203329, 134281, 203770,
                                                                       79856, 80136, 140896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 213178, 0, 3,
                                                                       203770, 134596, 204211,
                                                                       80136, 80416, 141316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 213766, 0, 3,
                                                                       204211, 134911, 204652,
                                                                       80416, 80696, 141736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 214354, 0, 3,
                                                                       204652, 135226, 205093,
                                                                       80696, 80976, 142156,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 214942, 0, 3,
                                                                       205534, 135856, 206122,
                                                                       81536, 81896, 142576,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 215698, 0, 3,
                                                                       206122, 136276, 206710,
                                                                       81896, 82256, 143116,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 216454, 0, 3,
                                                                       206710, 136696, 207298,
                                                                       82256, 82616, 143656,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 217210, 0, 3,
                                                                       207298, 137116, 207886,
                                                                       82616, 82976, 144196,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 217966, 0, 3,
                                                                       207886, 137536, 208474,
                                                                       82976, 83336, 144736,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 218722, 0, 3,
                                                                       208474, 137956, 209062,
                                                                       83336, 83696, 145276,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 219478, 0, 3,
                                                                       209062, 138376, 209650,
                                                                       83696, 84056, 145816,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 220234, 0, 3,
                                                                       210238, 139216, 210826,
                                                                       84776, 85136, 146356,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 220990, 0, 3,
                                                                       210826, 139636, 211414,
                                                                       85136, 85496, 146896,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 221746, 0, 3,
                                                                       211414, 140056, 212002,
                                                                       85496, 85856, 147436,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 222502, 0, 3,
                                                                       212002, 140476, 212590,
                                                                       85856, 86216, 147976,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 223258, 0, 3,
                                                                       212590, 140896, 213178,
                                                                       86216, 86576, 148516,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 224014, 0, 3,
                                                                       213178, 141316, 213766,
                                                                       86576, 86936, 149056,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 224770, 0, 3,
                                                                       213766, 141736, 214354,
                                                                       86936, 87296, 149596,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 225526, 0, 3,
                                                                       214942, 142576, 215698,
                                                                       88016, 88466, 150136,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 226471, 0, 3,
                                                                       215698, 143116, 216454,
                                                                       88466, 88916, 150811,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 227416, 0, 3,
                                                                       216454, 143656, 217210,
                                                                       88916, 89366, 151486,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 228361, 0, 3,
                                                                       217210, 144196, 217966,
                                                                       89366, 89816, 152161,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 229306, 0, 3,
                                                                       217966, 144736, 218722,
                                                                       89816, 90266, 152836,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 230251, 0, 3,
                                                                       218722, 145276, 219478,
                                                                       90266, 90716, 153511,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 231196, 0, 3,
                                                                       220234, 146356, 220990,
                                                                       91616, 92066, 154186,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 232141, 0, 3,
                                                                       220990, 146896, 221746,
                                                                       92066, 92516, 154861,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 233086, 0, 3,
                                                                       221746, 147436, 222502,
                                                                       92516, 92966, 155536,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 234031, 0, 3,
                                                                       222502, 147976, 223258,
                                                                       92966, 93416, 156211,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 234976, 0, 3,
                                                                       223258, 148516, 224014,
                                                                       93416, 93866, 156886,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 235921, 0, 3,
                                                                       224014, 149056, 224770,
                                                                       93866, 94316, 157561,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 236866, 0, 3,
                                                                       225526, 150136, 226471,
                                                                       95216, 95766, 158236,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 238021, 0, 3,
                                                                       226471, 150811, 227416,
                                                                       95766, 96316, 159061,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 239176, 0, 3,
                                                                       227416, 151486, 228361,
                                                                       96316, 96866, 159886,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 240331, 0, 3,
                                                                       228361, 152161, 229306,
                                                                       96866, 97416, 160711,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 241486, 0, 3,
                                                                       229306, 152836, 230251,
                                                                       97416, 97966, 161536,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 242641, 0, 3,
                                                                       231196, 154186, 232141,
                                                                       99066, 99616, 162361,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 243796, 0, 3,
                                                                       232141, 154861, 233086,
                                                                       99616, 100166, 163186,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 244951, 0, 3,
                                                                       233086, 155536, 234031,
                                                                       100166, 100716, 164011,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 246106, 0, 3,
                                                                       234031, 156211, 234976,
                                                                       100716, 101266, 164836,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 247261, 0, 3,
                                                                       234976, 156886, 235921,
                                                                       101266, 101816, 165661,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 248416, 0, 3,
                                                                       236866, 158236, 238021,
                                                                       102916, 103576, 166486,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 249802, 0, 3,
                                                                       238021, 159061, 239176,
                                                                       103576, 104236, 167476,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 251188, 0, 3,
                                                                       239176, 159886, 240331,
                                                                       104236, 104896, 168466,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 252574, 0, 3,
                                                                       240331, 160711, 241486,
                                                                       104896, 105556, 169456,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 253960, 0, 3,
                                                                       242641, 162361, 243796,
                                                                       106876, 107536, 170446,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 255346, 0, 3,
                                                                       243796, 163186, 244951,
                                                                       107536, 108196, 171436,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 256732, 0, 3,
                                                                       244951, 164011, 246106,
                                                                       108196, 108856, 172426,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 258118, 0, 3,
                                                                       246106, 164836, 247261,
                                                                       108856, 109516, 173416,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 259504, 0, 3,
                                                                       248416, 166486, 249802,
                                                                       110836, 111616, 174406,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 261142, 0, 3,
                                                                       249802, 167476, 251188,
                                                                       111616, 112396, 175576,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 262780, 0, 3,
                                                                       251188, 168466, 252574,
                                                                       112396, 113176, 176746,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 264418, 0, 3,
                                                                       253960, 170446, 255346,
                                                                       114736, 115516, 177916,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 266056, 0, 3,
                                                                       255346, 171436, 256732,
                                                                       115516, 116296, 179086,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 267694, 0, 3,
                                                                       256732, 172426, 258118,
                                                                       116296, 117076, 180256,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269332, 3, 118636,
                                                                       118651, 181468, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269360, 3, 118651,
                                                                       118666, 181489, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269388, 3, 118666,
                                                                       118681, 181510, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269416, 3, 118681,
                                                                       118696, 181531, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269444, 3, 118696,
                                                                       118711, 181552, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269472, 3, 118711,
                                                                       118726, 181573, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269500, 3, 118726,
                                                                       118741, 181594, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269528, 3, 118741,
                                                                       118756, 181615, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269556, 3, 118756,
                                                                       118771, 181636, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269584, 3, 118771,
                                                                       118786, 181657, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269612, 3, 118786,
                                                                       118801, 181678, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269640, 3, 118801,
                                                                       118816, 181699, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269668, 3, 118846,
                                                                       118861, 181762, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269696, 3, 118861,
                                                                       118876, 181783, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269724, 3, 118876,
                                                                       118891, 181804, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269752, 3, 118891,
                                                                       118906, 181825, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269780, 3, 118906,
                                                                       118921, 181846, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269808, 3, 118921,
                                                                       118936, 181867, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269836, 3, 118936,
                                                                       118951, 181888, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269864, 3, 118951,
                                                                       118966, 181909, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269892, 3, 118966,
                                                                       118981, 181930, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269920, 3, 118981,
                                                                       118996, 181951, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269948, 3, 118996,
                                                                       119011, 181972, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 269976, 3, 119011,
                                                                       119026, 181993, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270004, 0, 3,
                                                                       269332, 181468, 269360,
                                                                       119056, 119101, 182140,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270088, 0, 3,
                                                                       269360, 181489, 269388,
                                                                       119101, 119146, 182203,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270172, 0, 3,
                                                                       269388, 181510, 269416,
                                                                       119146, 119191, 182266,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270256, 0, 3,
                                                                       269416, 181531, 269444,
                                                                       119191, 119236, 182329,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270340, 0, 3,
                                                                       269444, 181552, 269472,
                                                                       119236, 119281, 182392,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270424, 0, 3,
                                                                       269472, 181573, 269500,
                                                                       119281, 119326, 182455,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270508, 0, 3,
                                                                       269500, 181594, 269528,
                                                                       119326, 119371, 182518,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270592, 0, 3,
                                                                       269528, 181615, 269556,
                                                                       119371, 119416, 182581,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270676, 0, 3,
                                                                       269556, 181636, 269584,
                                                                       119416, 119461, 182644,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270760, 0, 3,
                                                                       269584, 181657, 269612,
                                                                       119461, 119506, 182707,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270844, 0, 3,
                                                                       269612, 181678, 269640,
                                                                       119506, 119551, 182770,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 270928, 0, 3,
                                                                       269668, 181762, 269696,
                                                                       119641, 119686, 182959,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 271012, 0, 3,
                                                                       269696, 181783, 269724,
                                                                       119686, 119731, 183022,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 271096, 0, 3,
                                                                       269724, 181804, 269752,
                                                                       119731, 119776, 183085,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 271180, 0, 3,
                                                                       269752, 181825, 269780,
                                                                       119776, 119821, 183148,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 271264, 0, 3,
                                                                       269780, 181846, 269808,
                                                                       119821, 119866, 183211,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 271348, 0, 3,
                                                                       269808, 181867, 269836,
                                                                       119866, 119911, 183274,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 271432, 0, 3,
                                                                       269836, 181888, 269864,
                                                                       119911, 119956, 183337,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 271516, 0, 3,
                                                                       269864, 181909, 269892,
                                                                       119956, 120001, 183400,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 271600, 0, 3,
                                                                       269892, 181930, 269920,
                                                                       120001, 120046, 183463,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 271684, 0, 3,
                                                                       269920, 181951, 269948,
                                                                       120046, 120091, 183526,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 271768, 0, 3,
                                                                       269948, 181972, 269976,
                                                                       120091, 120136, 183589,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 271852, 0, 3,
                                                                       270004, 182140, 270088,
                                                                       120226, 120316, 183904,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 272020, 0, 3,
                                                                       270088, 182203, 270172,
                                                                       120316, 120406, 184030,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 272188, 0, 3,
                                                                       270172, 182266, 270256,
                                                                       120406, 120496, 184156,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 272356, 0, 3,
                                                                       270256, 182329, 270340,
                                                                       120496, 120586, 184282,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 272524, 0, 3,
                                                                       270340, 182392, 270424,
                                                                       120586, 120676, 184408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 272692, 0, 3,
                                                                       270424, 182455, 270508,
                                                                       120676, 120766, 184534,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 272860, 0, 3,
                                                                       270508, 182518, 270592,
                                                                       120766, 120856, 184660,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 273028, 0, 3,
                                                                       270592, 182581, 270676,
                                                                       120856, 120946, 184786,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 273196, 0, 3,
                                                                       270676, 182644, 270760,
                                                                       120946, 121036, 184912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 273364, 0, 3,
                                                                       270760, 182707, 270844,
                                                                       121036, 121126, 185038,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 273532, 0, 3,
                                                                       270928, 182959, 271012,
                                                                       121306, 121396, 185416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 273700, 0, 3,
                                                                       271012, 183022, 271096,
                                                                       121396, 121486, 185542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 273868, 0, 3,
                                                                       271096, 183085, 271180,
                                                                       121486, 121576, 185668,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 274036, 0, 3,
                                                                       271180, 183148, 271264,
                                                                       121576, 121666, 185794,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 274204, 0, 3,
                                                                       271264, 183211, 271348,
                                                                       121666, 121756, 185920,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 274372, 0, 3,
                                                                       271348, 183274, 271432,
                                                                       121756, 121846, 186046,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 274540, 0, 3,
                                                                       271432, 183337, 271516,
                                                                       121846, 121936, 186172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 274708, 0, 3,
                                                                       271516, 183400, 271600,
                                                                       121936, 122026, 186298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 274876, 0, 3,
                                                                       271600, 183463, 271684,
                                                                       122026, 122116, 186424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 275044, 0, 3,
                                                                       271684, 183526, 271768,
                                                                       122116, 122206, 186550,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 275212, 0, 3,
                                                                       271852, 183904, 272020,
                                                                       122386, 122536, 187096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 275492, 0, 3,
                                                                       272020, 184030, 272188,
                                                                       122536, 122686, 187306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 275772, 0, 3,
                                                                       272188, 184156, 272356,
                                                                       122686, 122836, 187516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 276052, 0, 3,
                                                                       272356, 184282, 272524,
                                                                       122836, 122986, 187726,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 276332, 0, 3,
                                                                       272524, 184408, 272692,
                                                                       122986, 123136, 187936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 276612, 0, 3,
                                                                       272692, 184534, 272860,
                                                                       123136, 123286, 188146,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 276892, 0, 3,
                                                                       272860, 184660, 273028,
                                                                       123286, 123436, 188356,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 277172, 0, 3,
                                                                       273028, 184786, 273196,
                                                                       123436, 123586, 188566,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 277452, 0, 3,
                                                                       273196, 184912, 273364,
                                                                       123586, 123736, 188776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 277732, 0, 3,
                                                                       273532, 185416, 273700,
                                                                       124036, 124186, 189406,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 278012, 0, 3,
                                                                       273700, 185542, 273868,
                                                                       124186, 124336, 189616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 278292, 0, 3,
                                                                       273868, 185668, 274036,
                                                                       124336, 124486, 189826,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 278572, 0, 3,
                                                                       274036, 185794, 274204,
                                                                       124486, 124636, 190036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 278852, 0, 3,
                                                                       274204, 185920, 274372,
                                                                       124636, 124786, 190246,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 279132, 0, 3,
                                                                       274372, 186046, 274540,
                                                                       124786, 124936, 190456,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 279412, 0, 3,
                                                                       274540, 186172, 274708,
                                                                       124936, 125086, 190666,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 279692, 0, 3,
                                                                       274708, 186298, 274876,
                                                                       125086, 125236, 190876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 279972, 0, 3,
                                                                       274876, 186424, 275044,
                                                                       125236, 125386, 191086,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 280252, 0, 3,
                                                                       275212, 187096, 275492,
                                                                       125686, 125911, 191926,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 280672, 0, 3,
                                                                       275492, 187306, 275772,
                                                                       125911, 126136, 192241,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 281092, 0, 3,
                                                                       275772, 187516, 276052,
                                                                       126136, 126361, 192556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 281512, 0, 3,
                                                                       276052, 187726, 276332,
                                                                       126361, 126586, 192871,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 281932, 0, 3,
                                                                       276332, 187936, 276612,
                                                                       126586, 126811, 193186,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 282352, 0, 3,
                                                                       276612, 188146, 276892,
                                                                       126811, 127036, 193501,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 282772, 0, 3,
                                                                       276892, 188356, 277172,
                                                                       127036, 127261, 193816,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 283192, 0, 3,
                                                                       277172, 188566, 277452,
                                                                       127261, 127486, 194131,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 283612, 0, 3,
                                                                       277732, 189406, 278012,
                                                                       127936, 128161, 195076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 284032, 0, 3,
                                                                       278012, 189616, 278292,
                                                                       128161, 128386, 195391,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 284452, 0, 3,
                                                                       278292, 189826, 278572,
                                                                       128386, 128611, 195706,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 284872, 0, 3,
                                                                       278572, 190036, 278852,
                                                                       128611, 128836, 196021,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 285292, 0, 3,
                                                                       278852, 190246, 279132,
                                                                       128836, 129061, 196336,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 285712, 0, 3,
                                                                       279132, 190456, 279412,
                                                                       129061, 129286, 196651,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 286132, 0, 3,
                                                                       279412, 190666, 279692,
                                                                       129286, 129511, 196966,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 286552, 0, 3,
                                                                       279692, 190876, 279972,
                                                                       129511, 129736, 197281,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 286972, 0, 3,
                                                                       280252, 191926, 280672,
                                                                       130186, 130501, 198478,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 287560, 0, 3,
                                                                       280672, 192241, 281092,
                                                                       130501, 130816, 198919,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 288148, 0, 3,
                                                                       281092, 192556, 281512,
                                                                       130816, 131131, 199360,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 288736, 0, 3,
                                                                       281512, 192871, 281932,
                                                                       131131, 131446, 199801,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 289324, 0, 3,
                                                                       281932, 193186, 282352,
                                                                       131446, 131761, 200242,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 289912, 0, 3,
                                                                       282352, 193501, 282772,
                                                                       131761, 132076, 200683,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 290500, 0, 3,
                                                                       282772, 193816, 283192,
                                                                       132076, 132391, 201124,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 291088, 0, 3,
                                                                       283612, 195076, 284032,
                                                                       133021, 133336, 202447,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 291676, 0, 3,
                                                                       284032, 195391, 284452,
                                                                       133336, 133651, 202888,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 292264, 0, 3,
                                                                       284452, 195706, 284872,
                                                                       133651, 133966, 203329,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 292852, 0, 3,
                                                                       284872, 196021, 285292,
                                                                       133966, 134281, 203770,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 293440, 0, 3,
                                                                       285292, 196336, 285712,
                                                                       134281, 134596, 204211,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 294028, 0, 3,
                                                                       285712, 196651, 286132,
                                                                       134596, 134911, 204652,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 294616, 0, 3,
                                                                       286132, 196966, 286552,
                                                                       134911, 135226, 205093,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 295204, 0, 3,
                                                                       286972, 198478, 287560,
                                                                       135856, 136276, 206710,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 295988, 0, 3,
                                                                       287560, 198919, 288148,
                                                                       136276, 136696, 207298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 296772, 0, 3,
                                                                       288148, 199360, 288736,
                                                                       136696, 137116, 207886,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 297556, 0, 3,
                                                                       288736, 199801, 289324,
                                                                       137116, 137536, 208474,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 298340, 0, 3,
                                                                       289324, 200242, 289912,
                                                                       137536, 137956, 209062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 299124, 0, 3,
                                                                       289912, 200683, 290500,
                                                                       137956, 138376, 209650,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 299908, 0, 3,
                                                                       291088, 202447, 291676,
                                                                       139216, 139636, 211414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 300692, 0, 3,
                                                                       291676, 202888, 292264,
                                                                       139636, 140056, 212002,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 301476, 0, 3,
                                                                       292264, 203329, 292852,
                                                                       140056, 140476, 212590,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 302260, 0, 3,
                                                                       292852, 203770, 293440,
                                                                       140476, 140896, 213178,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 303044, 0, 3,
                                                                       293440, 204211, 294028,
                                                                       140896, 141316, 213766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 303828, 0, 3,
                                                                       294028, 204652, 294616,
                                                                       141316, 141736, 214354,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 304612, 0, 3,
                                                                       295204, 206710, 295988,
                                                                       142576, 143116, 216454,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 305620, 0, 3,
                                                                       295988, 207298, 296772,
                                                                       143116, 143656, 217210,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 306628, 0, 3,
                                                                       296772, 207886, 297556,
                                                                       143656, 144196, 217966,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 307636, 0, 3,
                                                                       297556, 208474, 298340,
                                                                       144196, 144736, 218722,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 308644, 0, 3,
                                                                       298340, 209062, 299124,
                                                                       144736, 145276, 219478,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 309652, 0, 3,
                                                                       299908, 211414, 300692,
                                                                       146356, 146896, 221746,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 310660, 0, 3,
                                                                       300692, 212002, 301476,
                                                                       146896, 147436, 222502,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 311668, 0, 3,
                                                                       301476, 212590, 302260,
                                                                       147436, 147976, 223258,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 312676, 0, 3,
                                                                       302260, 213178, 303044,
                                                                       147976, 148516, 224014,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 313684, 0, 3,
                                                                       303044, 213766, 303828,
                                                                       148516, 149056, 224770,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 314692, 0, 3,
                                                                       304612, 216454, 305620,
                                                                       150136, 150811, 227416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 315952, 0, 3,
                                                                       305620, 217210, 306628,
                                                                       150811, 151486, 228361,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 317212, 0, 3,
                                                                       306628, 217966, 307636,
                                                                       151486, 152161, 229306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 318472, 0, 3,
                                                                       307636, 218722, 308644,
                                                                       152161, 152836, 230251,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 319732, 0, 3,
                                                                       309652, 221746, 310660,
                                                                       154186, 154861, 233086,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 320992, 0, 3,
                                                                       310660, 222502, 311668,
                                                                       154861, 155536, 234031,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 322252, 0, 3,
                                                                       311668, 223258, 312676,
                                                                       155536, 156211, 234976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 323512, 0, 3,
                                                                       312676, 224014, 313684,
                                                                       156211, 156886, 235921,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 324772, 0, 3,
                                                                       314692, 227416, 315952,
                                                                       158236, 159061, 239176,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 326312, 0, 3,
                                                                       315952, 228361, 317212,
                                                                       159061, 159886, 240331,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 327852, 0, 3,
                                                                       317212, 229306, 318472,
                                                                       159886, 160711, 241486,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 329392, 0, 3,
                                                                       319732, 233086, 320992,
                                                                       162361, 163186, 244951,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 330932, 0, 3,
                                                                       320992, 234031, 322252,
                                                                       163186, 164011, 246106,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 332472, 0, 3,
                                                                       322252, 234976, 323512,
                                                                       164011, 164836, 247261,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 334012, 0, 3,
                                                                       324772, 239176, 326312,
                                                                       166486, 167476, 251188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 335860, 0, 3,
                                                                       326312, 240331, 327852,
                                                                       167476, 168466, 252574,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 337708, 0, 3,
                                                                       329392, 244951, 330932,
                                                                       170446, 171436, 256732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 339556, 0, 3,
                                                                       330932, 246106, 332472,
                                                                       171436, 172426, 258118,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 341404, 0, 3,
                                                                       334012, 251188, 335860,
                                                                       174406, 175576, 262780,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 343588, 0, 3,
                                                                       337708, 256732, 339556,
                                                                       177916, 179086, 267694,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345772, 3, 181426,
                                                                       181447, 269332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345808, 3, 181447,
                                                                       181468, 269360, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345844, 3, 181468,
                                                                       181489, 269388, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345880, 3, 181489,
                                                                       181510, 269416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345916, 3, 181510,
                                                                       181531, 269444, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345952, 3, 181531,
                                                                       181552, 269472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 345988, 3, 181552,
                                                                       181573, 269500, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346024, 3, 181573,
                                                                       181594, 269528, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346060, 3, 181594,
                                                                       181615, 269556, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346096, 3, 181615,
                                                                       181636, 269584, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346132, 3, 181636,
                                                                       181657, 269612, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346168, 3, 181657,
                                                                       181678, 269640, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346204, 3, 181720,
                                                                       181741, 269668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346240, 3, 181741,
                                                                       181762, 269696, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346276, 3, 181762,
                                                                       181783, 269724, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346312, 3, 181783,
                                                                       181804, 269752, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346348, 3, 181804,
                                                                       181825, 269780, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346384, 3, 181825,
                                                                       181846, 269808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346420, 3, 181846,
                                                                       181867, 269836, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346456, 3, 181867,
                                                                       181888, 269864, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346492, 3, 181888,
                                                                       181909, 269892, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346528, 3, 181909,
                                                                       181930, 269920, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346564, 3, 181930,
                                                                       181951, 269948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 346600, 3, 181951,
                                                                       181972, 269976, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 346636, 0, 3,
                                                                       345772, 269332, 345808,
                                                                       182014, 182077, 270004,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 346744, 0, 3,
                                                                       345808, 269360, 345844,
                                                                       182077, 182140, 270088,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 346852, 0, 3,
                                                                       345844, 269388, 345880,
                                                                       182140, 182203, 270172,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 346960, 0, 3,
                                                                       345880, 269416, 345916,
                                                                       182203, 182266, 270256,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 347068, 0, 3,
                                                                       345916, 269444, 345952,
                                                                       182266, 182329, 270340,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 347176, 0, 3,
                                                                       345952, 269472, 345988,
                                                                       182329, 182392, 270424,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 347284, 0, 3,
                                                                       345988, 269500, 346024,
                                                                       182392, 182455, 270508,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 347392, 0, 3,
                                                                       346024, 269528, 346060,
                                                                       182455, 182518, 270592,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 347500, 0, 3,
                                                                       346060, 269556, 346096,
                                                                       182518, 182581, 270676,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 347608, 0, 3,
                                                                       346096, 269584, 346132,
                                                                       182581, 182644, 270760,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 347716, 0, 3,
                                                                       346132, 269612, 346168,
                                                                       182644, 182707, 270844,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 347824, 0, 3,
                                                                       346204, 269668, 346240,
                                                                       182833, 182896, 270928,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 347932, 0, 3,
                                                                       346240, 269696, 346276,
                                                                       182896, 182959, 271012,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 348040, 0, 3,
                                                                       346276, 269724, 346312,
                                                                       182959, 183022, 271096,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 348148, 0, 3,
                                                                       346312, 269752, 346348,
                                                                       183022, 183085, 271180,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 348256, 0, 3,
                                                                       346348, 269780, 346384,
                                                                       183085, 183148, 271264,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 348364, 0, 3,
                                                                       346384, 269808, 346420,
                                                                       183148, 183211, 271348,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 348472, 0, 3,
                                                                       346420, 269836, 346456,
                                                                       183211, 183274, 271432,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 348580, 0, 3,
                                                                       346456, 269864, 346492,
                                                                       183274, 183337, 271516,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 348688, 0, 3,
                                                                       346492, 269892, 346528,
                                                                       183337, 183400, 271600,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 348796, 0, 3,
                                                                       346528, 269920, 346564,
                                                                       183400, 183463, 271684,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 348904, 0, 3,
                                                                       346564, 269948, 346600,
                                                                       183463, 183526, 271768,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 349012, 0, 3,
                                                                       346636, 270004, 346744,
                                                                       183652, 183778, 271852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 349228, 0, 3,
                                                                       346744, 270088, 346852,
                                                                       183778, 183904, 272020,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 349444, 0, 3,
                                                                       346852, 270172, 346960,
                                                                       183904, 184030, 272188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 349660, 0, 3,
                                                                       346960, 270256, 347068,
                                                                       184030, 184156, 272356,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 349876, 0, 3,
                                                                       347068, 270340, 347176,
                                                                       184156, 184282, 272524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 350092, 0, 3,
                                                                       347176, 270424, 347284,
                                                                       184282, 184408, 272692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 350308, 0, 3,
                                                                       347284, 270508, 347392,
                                                                       184408, 184534, 272860,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 350524, 0, 3,
                                                                       347392, 270592, 347500,
                                                                       184534, 184660, 273028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 350740, 0, 3,
                                                                       347500, 270676, 347608,
                                                                       184660, 184786, 273196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 350956, 0, 3,
                                                                       347608, 270760, 347716,
                                                                       184786, 184912, 273364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 351172, 0, 3,
                                                                       347824, 270928, 347932,
                                                                       185164, 185290, 273532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 351388, 0, 3,
                                                                       347932, 271012, 348040,
                                                                       185290, 185416, 273700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 351604, 0, 3,
                                                                       348040, 271096, 348148,
                                                                       185416, 185542, 273868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 351820, 0, 3,
                                                                       348148, 271180, 348256,
                                                                       185542, 185668, 274036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 352036, 0, 3,
                                                                       348256, 271264, 348364,
                                                                       185668, 185794, 274204,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 352252, 0, 3,
                                                                       348364, 271348, 348472,
                                                                       185794, 185920, 274372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 352468, 0, 3,
                                                                       348472, 271432, 348580,
                                                                       185920, 186046, 274540,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 352684, 0, 3,
                                                                       348580, 271516, 348688,
                                                                       186046, 186172, 274708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 352900, 0, 3,
                                                                       348688, 271600, 348796,
                                                                       186172, 186298, 274876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 353116, 0, 3,
                                                                       348796, 271684, 348904,
                                                                       186298, 186424, 275044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 353332, 0, 3,
                                                                       349012, 271852, 349228,
                                                                       186676, 186886, 275212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 353692, 0, 3,
                                                                       349228, 272020, 349444,
                                                                       186886, 187096, 275492,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 354052, 0, 3,
                                                                       349444, 272188, 349660,
                                                                       187096, 187306, 275772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 354412, 0, 3,
                                                                       349660, 272356, 349876,
                                                                       187306, 187516, 276052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 354772, 0, 3,
                                                                       349876, 272524, 350092,
                                                                       187516, 187726, 276332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 355132, 0, 3,
                                                                       350092, 272692, 350308,
                                                                       187726, 187936, 276612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 355492, 0, 3,
                                                                       350308, 272860, 350524,
                                                                       187936, 188146, 276892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 355852, 0, 3,
                                                                       350524, 273028, 350740,
                                                                       188146, 188356, 277172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 356212, 0, 3,
                                                                       350740, 273196, 350956,
                                                                       188356, 188566, 277452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 356572, 0, 3,
                                                                       351172, 273532, 351388,
                                                                       188986, 189196, 277732,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 356932, 0, 3,
                                                                       351388, 273700, 351604,
                                                                       189196, 189406, 278012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 357292, 0, 3,
                                                                       351604, 273868, 351820,
                                                                       189406, 189616, 278292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 357652, 0, 3,
                                                                       351820, 274036, 352036,
                                                                       189616, 189826, 278572,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 358012, 0, 3,
                                                                       352036, 274204, 352252,
                                                                       189826, 190036, 278852,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 358372, 0, 3,
                                                                       352252, 274372, 352468,
                                                                       190036, 190246, 279132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 358732, 0, 3,
                                                                       352468, 274540, 352684,
                                                                       190246, 190456, 279412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 359092, 0, 3,
                                                                       352684, 274708, 352900,
                                                                       190456, 190666, 279692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 359452, 0, 3,
                                                                       352900, 274876, 353116,
                                                                       190666, 190876, 279972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 359812, 0, 3,
                                                                       353332, 275212, 353692,
                                                                       191296, 191611, 280252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 360352, 0, 3,
                                                                       353692, 275492, 354052,
                                                                       191611, 191926, 280672,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 360892, 0, 3,
                                                                       354052, 275772, 354412,
                                                                       191926, 192241, 281092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 361432, 0, 3,
                                                                       354412, 276052, 354772,
                                                                       192241, 192556, 281512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 361972, 0, 3,
                                                                       354772, 276332, 355132,
                                                                       192556, 192871, 281932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 362512, 0, 3,
                                                                       355132, 276612, 355492,
                                                                       192871, 193186, 282352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 363052, 0, 3,
                                                                       355492, 276892, 355852,
                                                                       193186, 193501, 282772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 363592, 0, 3,
                                                                       355852, 277172, 356212,
                                                                       193501, 193816, 283192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 364132, 0, 3,
                                                                       356572, 277732, 356932,
                                                                       194446, 194761, 283612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 364672, 0, 3,
                                                                       356932, 278012, 357292,
                                                                       194761, 195076, 284032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 365212, 0, 3,
                                                                       357292, 278292, 357652,
                                                                       195076, 195391, 284452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 365752, 0, 3,
                                                                       357652, 278572, 358012,
                                                                       195391, 195706, 284872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 366292, 0, 3,
                                                                       358012, 278852, 358372,
                                                                       195706, 196021, 285292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 366832, 0, 3,
                                                                       358372, 279132, 358732,
                                                                       196021, 196336, 285712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 367372, 0, 3,
                                                                       358732, 279412, 359092,
                                                                       196336, 196651, 286132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 367912, 0, 3,
                                                                       359092, 279692, 359452,
                                                                       196651, 196966, 286552,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 368452, 0, 3,
                                                                       359812, 280252, 360352,
                                                                       197596, 198037, 286972,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 369208, 0, 3,
                                                                       360352, 280672, 360892,
                                                                       198037, 198478, 287560,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 369964, 0, 3,
                                                                       360892, 281092, 361432,
                                                                       198478, 198919, 288148,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 370720, 0, 3,
                                                                       361432, 281512, 361972,
                                                                       198919, 199360, 288736,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 371476, 0, 3,
                                                                       361972, 281932, 362512,
                                                                       199360, 199801, 289324,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 372232, 0, 3,
                                                                       362512, 282352, 363052,
                                                                       199801, 200242, 289912,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 372988, 0, 3,
                                                                       363052, 282772, 363592,
                                                                       200242, 200683, 290500,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 373744, 0, 3,
                                                                       364132, 283612, 364672,
                                                                       201565, 202006, 291088,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 374500, 0, 3,
                                                                       364672, 284032, 365212,
                                                                       202006, 202447, 291676,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 375256, 0, 3,
                                                                       365212, 284452, 365752,
                                                                       202447, 202888, 292264,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 376012, 0, 3,
                                                                       365752, 284872, 366292,
                                                                       202888, 203329, 292852,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 376768, 0, 3,
                                                                       366292, 285292, 366832,
                                                                       203329, 203770, 293440,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 377524, 0, 3,
                                                                       366832, 285712, 367372,
                                                                       203770, 204211, 294028,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 378280, 0, 3,
                                                                       367372, 286132, 367912,
                                                                       204211, 204652, 294616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 379036, 0, 3,
                                                                       368452, 286972, 369208,
                                                                       205534, 206122, 295204,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 380044, 0, 3,
                                                                       369208, 287560, 369964,
                                                                       206122, 206710, 295988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 381052, 0, 3,
                                                                       369964, 288148, 370720,
                                                                       206710, 207298, 296772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 382060, 0, 3,
                                                                       370720, 288736, 371476,
                                                                       207298, 207886, 297556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 383068, 0, 3,
                                                                       371476, 289324, 372232,
                                                                       207886, 208474, 298340,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 384076, 0, 3,
                                                                       372232, 289912, 372988,
                                                                       208474, 209062, 299124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 385084, 0, 3,
                                                                       373744, 291088, 374500,
                                                                       210238, 210826, 299908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 386092, 0, 3,
                                                                       374500, 291676, 375256,
                                                                       210826, 211414, 300692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 387100, 0, 3,
                                                                       375256, 292264, 376012,
                                                                       211414, 212002, 301476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 388108, 0, 3,
                                                                       376012, 292852, 376768,
                                                                       212002, 212590, 302260,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 389116, 0, 3,
                                                                       376768, 293440, 377524,
                                                                       212590, 213178, 303044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 390124, 0, 3,
                                                                       377524, 294028, 378280,
                                                                       213178, 213766, 303828,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 391132, 0, 3,
                                                                       379036, 295204, 380044,
                                                                       214942, 215698, 304612,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 392428, 0, 3,
                                                                       380044, 295988, 381052,
                                                                       215698, 216454, 305620,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 393724, 0, 3,
                                                                       381052, 296772, 382060,
                                                                       216454, 217210, 306628,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 395020, 0, 3,
                                                                       382060, 297556, 383068,
                                                                       217210, 217966, 307636,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 396316, 0, 3,
                                                                       383068, 298340, 384076,
                                                                       217966, 218722, 308644,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 397612, 0, 3,
                                                                       385084, 299908, 386092,
                                                                       220234, 220990, 309652,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 398908, 0, 3,
                                                                       386092, 300692, 387100,
                                                                       220990, 221746, 310660,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 400204, 0, 3,
                                                                       387100, 301476, 388108,
                                                                       221746, 222502, 311668,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 401500, 0, 3,
                                                                       388108, 302260, 389116,
                                                                       222502, 223258, 312676,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 402796, 0, 3,
                                                                       389116, 303044, 390124,
                                                                       223258, 224014, 313684,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 404092, 0, 3,
                                                                       391132, 304612, 392428,
                                                                       225526, 226471, 314692,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 405712, 0, 3,
                                                                       392428, 305620, 393724,
                                                                       226471, 227416, 315952,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 407332, 0, 3,
                                                                       393724, 306628, 395020,
                                                                       227416, 228361, 317212,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 408952, 0, 3,
                                                                       395020, 307636, 396316,
                                                                       228361, 229306, 318472,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 410572, 0, 3,
                                                                       397612, 309652, 398908,
                                                                       231196, 232141, 319732,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 412192, 0, 3,
                                                                       398908, 310660, 400204,
                                                                       232141, 233086, 320992,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 413812, 0, 3,
                                                                       400204, 311668, 401500,
                                                                       233086, 234031, 322252,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 415432, 0, 3,
                                                                       401500, 312676, 402796,
                                                                       234031, 234976, 323512,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 417052, 0, 3,
                                                                       404092, 314692, 405712,
                                                                       236866, 238021, 324772,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 419032, 0, 3,
                                                                       405712, 315952, 407332,
                                                                       238021, 239176, 326312,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 421012, 0, 3,
                                                                       407332, 317212, 408952,
                                                                       239176, 240331, 327852,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 422992, 0, 3,
                                                                       410572, 319732, 412192,
                                                                       242641, 243796, 329392,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 424972, 0, 3,
                                                                       412192, 320992, 413812,
                                                                       243796, 244951, 330932,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 426952, 0, 3,
                                                                       413812, 322252, 415432,
                                                                       244951, 246106, 332472,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 428932, 0, 3,
                                                                       417052, 324772, 419032,
                                                                       248416, 249802, 334012,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 431308, 0, 3,
                                                                       419032, 326312, 421012,
                                                                       249802, 251188, 335860,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 433684, 0, 3,
                                                                       422992, 329392, 424972,
                                                                       253960, 255346, 337708,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 436060, 0, 3,
                                                                       424972, 330932, 426952,
                                                                       255346, 256732, 339556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sok_three_center_electron_repulsion_0(buffer, 438436, 0, 3,
                                                                       428932, 334012, 431308,
                                                                       259504, 261142, 341404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sok_three_center_electron_repulsion_0(buffer, 441244, 0, 3,
                                                                       433684, 337708, 436060,
                                                                       264418, 266056, 343588,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 444052, 379036, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 445480, 385084, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 446908, 391132, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 448744, 397612, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 450580, 404092, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 452875, 410572, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 455170, 417052, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 457975, 422992, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 460780, 428932, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 464146, 433684, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 467512, 438436, 2808, ncols);

                    simdfunc::contract_primitives(buffer, 471490, 441244, 2808, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 445060, 444052, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 446488, 445480, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 448204, 446908, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 450040, 448744, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 452200, 450580, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 454495, 452875, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 457150, 455170, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 459955, 457975, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 463156, 460780, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 466522, 464146, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 470320, 467512, 78, 1, nmax);

        simdtrf::transform_k_inner(buffer, 474298, 471490, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 475468, 445060, 448204, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 476728, 446488, 450040, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 477988, 448204, 452200, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 479608, 450040, 454495, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 481228, 452200, 457150, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 483253, 454495, 459955, 15, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 485278, 457150, 463156, 15, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 487753, 459955, 466522, 15, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 490228, 463156, 470320, 15, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 493198, 466522, 474298, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 496168, 475468, 477988, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 498688, 476728, 479608, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 501208, 477988, 481228, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 504448, 479608, 483253, 15, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 507688, 481228, 485278, 15, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 511738, 483253, 487753, 15, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 515788, 485278, 490228, 15, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 520738, 487753, 493198, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 525688, 496168, 501208, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 529888, 498688, 504448, 15, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 534088, 501208, 507688, 15, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 539488, 504448, 511738, 15, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 544888, 507688, 515788, 15, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 551638, 511738, 520738, 15, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 558388, 525688, 534088, 15, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 564688, 529888, 539488, 15, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 570988, 534088, 544888, 15, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 579088, 539488, 551638, 15, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 587188, 558388, 570988, 15, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 596008, 564688, 579088, 15, nmax);

        simdtrf::transform_i_inner(buffer, 604828, 596008, 21, 15, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 604828, 195, nmax);

        simdtrf::transform_i_inner(buffer, 604828, 587188, 21, 15, nmax);

        simdtrf::transform_h_outer(values + 2145 * nvalues + n * npairs, nvalues, buffer, 604828,
                                   195, nmax);
    }

    for (size_t m = 0; m < 4290; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
