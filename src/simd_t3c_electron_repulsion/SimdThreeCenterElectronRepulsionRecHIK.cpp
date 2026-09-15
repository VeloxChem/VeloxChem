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


#include "SimdThreeCenterElectronRepulsionRecHIK.hpp"

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
compute_hik_three_center_electron_repulsion(double               *values,
                                            const size_t          npairs,
                                            const size_t          natoms,
                                            const CBasisFunction &a_function,
                                            const CBasisFunction &b_function,
                                            const CBasisFunction &c_function,
                                            const CSimdMatrix    &coordinates,
                                            const CSimdMatrix    &c_coordinates,
                                            CSimdMatrix          &buffer,
                                            const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_hik_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 306513, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2145 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 306513, 222030, 14538, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                simdfunc::compute_pair_exponent(buffer, coordinates, 6, nmax, mu);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 7, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17, 18}, ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 77, 0, 3, 8, 9,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 83, 0, 3, 9, 10,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 89, 0, 3, 10, 11,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 95, 0, 3, 11, 12,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 101, 0, 3, 12, 13,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 107, 0, 3, 13, 14,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 113, 0, 3, 14, 15,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 119, 0, 3, 15, 16,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 125, 0, 3, 16, 17,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 131, 0, 3, 17, 18,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 137, 0, 3, 18, 19,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 143, 0, 3, 19, 20,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 149, 0, 3, 20, 21,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 155, 0, 3, 21, 22,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 161, 0, 3, 22, 23,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 167, 0, 3, 23, 24,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 26, 29,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 29, 32,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 32, 35,
                                                                       89, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 35, 38,
                                                                       95, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 38, 41,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 41, 44,
                                                                       107, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 233, 0, 3, 44, 47,
                                                                       113, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 47, 50,
                                                                       119, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 50, 53,
                                                                       125, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 263, 0, 3, 53, 56,
                                                                       131, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 56, 59,
                                                                       137, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 283, 0, 3, 59, 62,
                                                                       143, 149, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 293, 0, 3, 62, 65,
                                                                       149, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 303, 0, 3, 65, 68,
                                                                       155, 161, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 313, 0, 3, 68, 71,
                                                                       161, 167, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 323, 0, 3, 77, 83,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 338, 0, 3, 83, 89,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 353, 0, 3, 89, 95,
                                                                       193, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 368, 0, 3, 95,
                                                                       101, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 383, 0, 3, 101,
                                                                       107, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 398, 0, 3, 107,
                                                                       113, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 413, 0, 3, 113,
                                                                       119, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 428, 0, 3, 119,
                                                                       125, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 443, 0, 3, 125,
                                                                       131, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 458, 0, 3, 131,
                                                                       137, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 473, 0, 3, 137,
                                                                       143, 273, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 488, 0, 3, 143,
                                                                       149, 283, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 503, 0, 3, 149,
                                                                       155, 293, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 518, 0, 3, 155,
                                                                       161, 303, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 533, 0, 3, 173,
                                                                       183, 323, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 554, 0, 3, 183,
                                                                       193, 338, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 575, 0, 3, 193,
                                                                       203, 353, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 596, 0, 3, 203,
                                                                       213, 368, 383, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 617, 0, 3, 213,
                                                                       223, 383, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 638, 0, 3, 223,
                                                                       233, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 659, 0, 3, 233,
                                                                       243, 413, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 680, 0, 3, 243,
                                                                       253, 428, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 701, 0, 3, 253,
                                                                       263, 443, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 722, 0, 3, 263,
                                                                       273, 458, 473, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 743, 0, 3, 273,
                                                                       283, 473, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 764, 0, 3, 283,
                                                                       293, 488, 503, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 785, 0, 3, 293,
                                                                       303, 503, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 806, 0, 3, 323,
                                                                       338, 533, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 834, 0, 3, 338,
                                                                       353, 554, 575, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 862, 0, 3, 353,
                                                                       368, 575, 596, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 890, 0, 3, 368,
                                                                       383, 596, 617, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 918, 0, 3, 383,
                                                                       398, 617, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 946, 0, 3, 398,
                                                                       413, 638, 659, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 974, 0, 3, 413,
                                                                       428, 659, 680, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 428,
                                                                       443, 680, 701, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 443,
                                                                       458, 701, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 458,
                                                                       473, 722, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 473,
                                                                       488, 743, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 488,
                                                                       503, 764, 785, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 533,
                                                                       554, 806, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1178, 0, 3, 554,
                                                                       575, 834, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1214, 0, 3, 575,
                                                                       596, 862, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1250, 0, 3, 596,
                                                                       617, 890, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1286, 0, 3, 617,
                                                                       638, 918, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1322, 0, 3, 638,
                                                                       659, 946, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1358, 0, 3, 659,
                                                                       680, 974, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1394, 0, 3, 680,
                                                                       701, 1002, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1430, 0, 3, 701,
                                                                       722, 1030, 1058, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1466, 0, 3, 722,
                                                                       743, 1058, 1086, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1502, 0, 3, 743,
                                                                       764, 1086, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1538, 0, 3, 806,
                                                                       834, 1142, 1178, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1583, 0, 3, 834,
                                                                       862, 1178, 1214, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1628, 0, 3, 862,
                                                                       890, 1214, 1250, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1673, 0, 3, 890,
                                                                       918, 1250, 1286, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1718, 0, 3, 918,
                                                                       946, 1286, 1322, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1763, 0, 3, 946,
                                                                       974, 1322, 1358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1808, 0, 3, 974,
                                                                       1002, 1358, 1394, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1853, 0, 3, 1002,
                                                                       1030, 1394, 1430, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1898, 0, 3, 1030,
                                                                       1058, 1430, 1466, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1943, 0, 3, 1058,
                                                                       1086, 1466, 1502, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1142,
                                                                       1178, 1538, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2043, 0, 3, 1178,
                                                                       1214, 1583, 1628, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2098, 0, 3, 1214,
                                                                       1250, 1628, 1673, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1250,
                                                                       1286, 1673, 1718, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2208, 0, 3, 1286,
                                                                       1322, 1718, 1763, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2263, 0, 3, 1322,
                                                                       1358, 1763, 1808, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2318, 0, 3, 1358,
                                                                       1394, 1808, 1853, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2373, 0, 3, 1394,
                                                                       1430, 1853, 1898, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2428, 0, 3, 1430,
                                                                       1466, 1898, 1943, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2483, 0, 3, 1538,
                                                                       1583, 1988, 2043, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2549, 0, 3, 1583,
                                                                       1628, 2043, 2098, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2615, 0, 3, 1628,
                                                                       1673, 2098, 2153, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2681, 0, 3, 1673,
                                                                       1718, 2153, 2208, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2747, 0, 3, 1718,
                                                                       1763, 2208, 2263, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2813, 0, 3, 1763,
                                                                       1808, 2263, 2318, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2879, 0, 3, 1808,
                                                                       1853, 2318, 2373, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2945, 0, 3, 1853,
                                                                       1898, 2373, 2428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3011, 0, 3, 1988,
                                                                       2043, 2483, 2549, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3089, 0, 3, 2043,
                                                                       2098, 2549, 2615, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3167, 0, 3, 2098,
                                                                       2153, 2615, 2681, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3245, 0, 3, 2153,
                                                                       2208, 2681, 2747, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3323, 0, 3, 2208,
                                                                       2263, 2747, 2813, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3401, 0, 3, 2263,
                                                                       2318, 2813, 2879, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3479, 0, 3, 2318,
                                                                       2373, 2879, 2945, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3557, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3560, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3563, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3566, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3569, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3572, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3575, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3578, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3581, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3584, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3587, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3590, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3593, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3596, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3599, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3602, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3605, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3608, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3611, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3620, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3629, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3638, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3647, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3656, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3665, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3674, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3683, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3692, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3701, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3710, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3719, 3, 22, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3728, 3, 23, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3737, 3, 24, 74,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3746, 3, 26, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3764, 3, 29, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3782, 3, 32, 89,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3800, 3, 35, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3818, 3, 38, 101,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3836, 3, 41, 107,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3854, 3, 44, 113,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3872, 3, 47, 119,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3890, 3, 50, 125,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3908, 3, 53, 131,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3926, 3, 56, 137,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3944, 3, 59, 143,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3962, 3, 62, 149,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3980, 3, 65, 155,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3998, 3, 68, 161,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4016, 3, 71, 167,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4034, 3, 77, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4064, 3, 83, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4094, 3, 89, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4124, 3, 95, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4154, 3, 101, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4184, 3, 107, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4214, 3, 113, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4244, 3, 119, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4274, 3, 125, 253,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4304, 3, 131, 263,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4334, 3, 137, 273,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4364, 3, 143, 283,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4394, 3, 149, 293,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4424, 3, 155, 303,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4454, 3, 161, 313,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4484, 3, 173, 323,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4529, 3, 183, 338,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4574, 3, 193, 353,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4619, 3, 203, 368,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4664, 3, 213, 383,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4709, 3, 223, 398,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4754, 3, 233, 413,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4799, 3, 243, 428,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4844, 3, 253, 443,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4889, 3, 263, 458,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4934, 3, 273, 473,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4979, 3, 283, 488,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5024, 3, 293, 503,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5069, 3, 303, 518,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5114, 3, 323, 533,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5177, 3, 338, 554,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5240, 3, 353, 575,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5303, 3, 368, 596,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5366, 3, 383, 617,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5429, 3, 398, 638,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5492, 3, 413, 659,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5555, 3, 428, 680,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5618, 3, 443, 701,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5681, 3, 458, 722,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5744, 3, 473, 743,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5807, 3, 488, 764,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5870, 3, 503, 785,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5933, 3, 533, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6017, 3, 554, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6101, 3, 575, 862,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6185, 3, 596, 890,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6269, 3, 617, 918,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6353, 3, 638, 946,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6437, 3, 659, 974,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6521, 3, 680,
                                                                       1002, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6605, 3, 701,
                                                                       1030, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6689, 3, 722,
                                                                       1058, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6773, 3, 743,
                                                                       1086, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6857, 3, 764,
                                                                       1114, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6941, 3, 806,
                                                                       1142, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7049, 3, 834,
                                                                       1178, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7157, 3, 862,
                                                                       1214, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7265, 3, 890,
                                                                       1250, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7373, 3, 918,
                                                                       1286, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7481, 3, 946,
                                                                       1322, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7589, 3, 974,
                                                                       1358, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7697, 3, 1002,
                                                                       1394, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7805, 3, 1030,
                                                                       1430, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7913, 3, 1058,
                                                                       1466, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8021, 3, 1086,
                                                                       1502, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8129, 3, 1142,
                                                                       1538, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8264, 3, 1178,
                                                                       1583, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8399, 3, 1214,
                                                                       1628, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8534, 3, 1250,
                                                                       1673, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8669, 3, 1286,
                                                                       1718, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8804, 3, 1322,
                                                                       1763, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8939, 3, 1358,
                                                                       1808, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9074, 3, 1394,
                                                                       1853, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9209, 3, 1430,
                                                                       1898, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9344, 3, 1466,
                                                                       1943, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9479, 3, 1538,
                                                                       1988, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9644, 3, 1583,
                                                                       2043, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9809, 3, 1628,
                                                                       2098, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9974, 3, 1673,
                                                                       2153, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10139, 3, 1718,
                                                                       2208, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10304, 3, 1763,
                                                                       2263, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10469, 3, 1808,
                                                                       2318, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10634, 3, 1853,
                                                                       2373, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10799, 3, 1898,
                                                                       2428, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10964, 3, 1988,
                                                                       2483, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11162, 3, 2043,
                                                                       2549, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11360, 3, 2098,
                                                                       2615, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11558, 3, 2153,
                                                                       2681, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11756, 3, 2208,
                                                                       2747, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11954, 3, 2263,
                                                                       2813, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12152, 3, 2318,
                                                                       2879, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12350, 3, 2373,
                                                                       2945, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 12548, 3, 2483,
                                                                       3011, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 12782, 3, 2549,
                                                                       3089, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13016, 3, 2615,
                                                                       3167, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13250, 3, 2681,
                                                                       3245, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13484, 3, 2747,
                                                                       3323, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13718, 3, 2813,
                                                                       3401, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13952, 3, 2879,
                                                                       3479, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14186, 3, 8, 9,
                                                                       3563, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14192, 3, 9, 10,
                                                                       3566, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14198, 3, 10, 11,
                                                                       3569, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14204, 3, 11, 12,
                                                                       3572, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14210, 3, 12, 13,
                                                                       3575, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14216, 3, 13, 14,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14222, 3, 14, 15,
                                                                       3581, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14228, 3, 15, 16,
                                                                       3584, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14234, 3, 16, 17,
                                                                       3587, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14240, 3, 17, 18,
                                                                       3590, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14246, 3, 18, 19,
                                                                       3593, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14252, 3, 19, 20,
                                                                       3596, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14258, 3, 20, 21,
                                                                       3599, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14264, 3, 21, 22,
                                                                       3602, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14270, 3, 22, 23,
                                                                       3605, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14276, 3, 23, 24,
                                                                       3608, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14282, 0, 3,
                                                                       14186, 3563, 14192, 3611,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14300, 0, 3,
                                                                       14192, 3566, 14198, 3620,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14318, 0, 3,
                                                                       14198, 3569, 14204, 3629,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14336, 0, 3,
                                                                       14204, 3572, 14210, 3638,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14354, 0, 3,
                                                                       14210, 3575, 14216, 3647,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14372, 0, 3,
                                                                       14216, 3578, 14222, 3656,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14390, 0, 3,
                                                                       14222, 3581, 14228, 3665,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14408, 0, 3,
                                                                       14228, 3584, 14234, 3674,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14426, 0, 3,
                                                                       14234, 3587, 14240, 3683,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14444, 0, 3,
                                                                       14240, 3590, 14246, 3692,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14462, 0, 3,
                                                                       14246, 3593, 14252, 3701,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14480, 0, 3,
                                                                       14252, 3596, 14258, 3710,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14498, 0, 3,
                                                                       14258, 3599, 14264, 3719,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14516, 0, 3,
                                                                       14264, 3602, 14270, 3728,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14534, 0, 3,
                                                                       14270, 3605, 14276, 3737,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14552, 0, 3,
                                                                       14282, 3611, 14300, 77,
                                                                       83, 3782, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14588, 0, 3,
                                                                       14300, 3620, 14318, 83,
                                                                       89, 3800, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14624, 0, 3,
                                                                       14318, 3629, 14336, 89,
                                                                       95, 3818, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14660, 0, 3,
                                                                       14336, 3638, 14354, 95,
                                                                       101, 3836, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14696, 0, 3,
                                                                       14354, 3647, 14372, 101,
                                                                       107, 3854, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14732, 0, 3,
                                                                       14372, 3656, 14390, 107,
                                                                       113, 3872, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14768, 0, 3,
                                                                       14390, 3665, 14408, 113,
                                                                       119, 3890, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14804, 0, 3,
                                                                       14408, 3674, 14426, 119,
                                                                       125, 3908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14840, 0, 3,
                                                                       14426, 3683, 14444, 125,
                                                                       131, 3926, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14876, 0, 3,
                                                                       14444, 3692, 14462, 131,
                                                                       137, 3944, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14912, 0, 3,
                                                                       14462, 3701, 14480, 137,
                                                                       143, 3962, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14948, 0, 3,
                                                                       14480, 3710, 14498, 143,
                                                                       149, 3980, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14984, 0, 3,
                                                                       14498, 3719, 14516, 149,
                                                                       155, 3998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15020, 0, 3,
                                                                       14516, 3728, 14534, 155,
                                                                       161, 4016, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15056, 0, 3,
                                                                       14552, 3782, 14588, 173,
                                                                       183, 4094, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15116, 0, 3,
                                                                       14588, 3800, 14624, 183,
                                                                       193, 4124, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15176, 0, 3,
                                                                       14624, 3818, 14660, 193,
                                                                       203, 4154, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15236, 0, 3,
                                                                       14660, 3836, 14696, 203,
                                                                       213, 4184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15296, 0, 3,
                                                                       14696, 3854, 14732, 213,
                                                                       223, 4214, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15356, 0, 3,
                                                                       14732, 3872, 14768, 223,
                                                                       233, 4244, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15416, 0, 3,
                                                                       14768, 3890, 14804, 233,
                                                                       243, 4274, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15476, 0, 3,
                                                                       14804, 3908, 14840, 243,
                                                                       253, 4304, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15536, 0, 3,
                                                                       14840, 3926, 14876, 253,
                                                                       263, 4334, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15596, 0, 3,
                                                                       14876, 3944, 14912, 263,
                                                                       273, 4364, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15656, 0, 3,
                                                                       14912, 3962, 14948, 273,
                                                                       283, 4394, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15716, 0, 3,
                                                                       14948, 3980, 14984, 283,
                                                                       293, 4424, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15776, 0, 3,
                                                                       14984, 3998, 15020, 293,
                                                                       303, 4454, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15836, 0, 3,
                                                                       15056, 4094, 15116, 323,
                                                                       338, 4574, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15926, 0, 3,
                                                                       15116, 4124, 15176, 338,
                                                                       353, 4619, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16016, 0, 3,
                                                                       15176, 4154, 15236, 353,
                                                                       368, 4664, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16106, 0, 3,
                                                                       15236, 4184, 15296, 368,
                                                                       383, 4709, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16196, 0, 3,
                                                                       15296, 4214, 15356, 383,
                                                                       398, 4754, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16286, 0, 3,
                                                                       15356, 4244, 15416, 398,
                                                                       413, 4799, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16376, 0, 3,
                                                                       15416, 4274, 15476, 413,
                                                                       428, 4844, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16466, 0, 3,
                                                                       15476, 4304, 15536, 428,
                                                                       443, 4889, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16556, 0, 3,
                                                                       15536, 4334, 15596, 443,
                                                                       458, 4934, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16646, 0, 3,
                                                                       15596, 4364, 15656, 458,
                                                                       473, 4979, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16736, 0, 3,
                                                                       15656, 4394, 15716, 473,
                                                                       488, 5024, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16826, 0, 3,
                                                                       15716, 4424, 15776, 488,
                                                                       503, 5069, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16916, 0, 3,
                                                                       15836, 4574, 15926, 533,
                                                                       554, 5240, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17042, 0, 3,
                                                                       15926, 4619, 16016, 554,
                                                                       575, 5303, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17168, 0, 3,
                                                                       16016, 4664, 16106, 575,
                                                                       596, 5366, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17294, 0, 3,
                                                                       16106, 4709, 16196, 596,
                                                                       617, 5429, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17420, 0, 3,
                                                                       16196, 4754, 16286, 617,
                                                                       638, 5492, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17546, 0, 3,
                                                                       16286, 4799, 16376, 638,
                                                                       659, 5555, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17672, 0, 3,
                                                                       16376, 4844, 16466, 659,
                                                                       680, 5618, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17798, 0, 3,
                                                                       16466, 4889, 16556, 680,
                                                                       701, 5681, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17924, 0, 3,
                                                                       16556, 4934, 16646, 701,
                                                                       722, 5744, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18050, 0, 3,
                                                                       16646, 4979, 16736, 722,
                                                                       743, 5807, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18176, 0, 3,
                                                                       16736, 5024, 16826, 743,
                                                                       764, 5870, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18302, 0, 3,
                                                                       16916, 5240, 17042, 806,
                                                                       834, 6101, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18470, 0, 3,
                                                                       17042, 5303, 17168, 834,
                                                                       862, 6185, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18638, 0, 3,
                                                                       17168, 5366, 17294, 862,
                                                                       890, 6269, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18806, 0, 3,
                                                                       17294, 5429, 17420, 890,
                                                                       918, 6353, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18974, 0, 3,
                                                                       17420, 5492, 17546, 918,
                                                                       946, 6437, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19142, 0, 3,
                                                                       17546, 5555, 17672, 946,
                                                                       974, 6521, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19310, 0, 3,
                                                                       17672, 5618, 17798, 974,
                                                                       1002, 6605, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19478, 0, 3,
                                                                       17798, 5681, 17924, 1002,
                                                                       1030, 6689, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19646, 0, 3,
                                                                       17924, 5744, 18050, 1030,
                                                                       1058, 6773, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19814, 0, 3,
                                                                       18050, 5807, 18176, 1058,
                                                                       1086, 6857, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19982, 0, 3,
                                                                       18302, 6101, 18470, 1142,
                                                                       1178, 7157, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20198, 0, 3,
                                                                       18470, 6185, 18638, 1178,
                                                                       1214, 7265, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20414, 0, 3,
                                                                       18638, 6269, 18806, 1214,
                                                                       1250, 7373, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20630, 0, 3,
                                                                       18806, 6353, 18974, 1250,
                                                                       1286, 7481, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20846, 0, 3,
                                                                       18974, 6437, 19142, 1286,
                                                                       1322, 7589, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21062, 0, 3,
                                                                       19142, 6521, 19310, 1322,
                                                                       1358, 7697, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21278, 0, 3,
                                                                       19310, 6605, 19478, 1358,
                                                                       1394, 7805, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21494, 0, 3,
                                                                       19478, 6689, 19646, 1394,
                                                                       1430, 7913, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21710, 0, 3,
                                                                       19646, 6773, 19814, 1430,
                                                                       1466, 8021, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21926, 0, 3,
                                                                       19982, 7157, 20198, 1538,
                                                                       1583, 8399, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 22196, 0, 3,
                                                                       20198, 7265, 20414, 1583,
                                                                       1628, 8534, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 22466, 0, 3,
                                                                       20414, 7373, 20630, 1628,
                                                                       1673, 8669, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 22736, 0, 3,
                                                                       20630, 7481, 20846, 1673,
                                                                       1718, 8804, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23006, 0, 3,
                                                                       20846, 7589, 21062, 1718,
                                                                       1763, 8939, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23276, 0, 3,
                                                                       21062, 7697, 21278, 1763,
                                                                       1808, 9074, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23546, 0, 3,
                                                                       21278, 7805, 21494, 1808,
                                                                       1853, 9209, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23816, 0, 3,
                                                                       21494, 7913, 21710, 1853,
                                                                       1898, 9344, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 24086, 0, 3,
                                                                       21926, 8399, 22196, 1988,
                                                                       2043, 9809, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 24416, 0, 3,
                                                                       22196, 8534, 22466, 2043,
                                                                       2098, 9974, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 24746, 0, 3,
                                                                       22466, 8669, 22736, 2098,
                                                                       2153, 10139, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 25076, 0, 3,
                                                                       22736, 8804, 23006, 2153,
                                                                       2208, 10304, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 25406, 0, 3,
                                                                       23006, 8939, 23276, 2208,
                                                                       2263, 10469, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 25736, 0, 3,
                                                                       23276, 9074, 23546, 2263,
                                                                       2318, 10634, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 26066, 0, 3,
                                                                       23546, 9209, 23816, 2318,
                                                                       2373, 10799, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 26396, 0, 3,
                                                                       24086, 9809, 24416, 2483,
                                                                       2549, 11360, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 26792, 0, 3,
                                                                       24416, 9974, 24746, 2549,
                                                                       2615, 11558, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 27188, 0, 3,
                                                                       24746, 10139, 25076, 2615,
                                                                       2681, 11756, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 27584, 0, 3,
                                                                       25076, 10304, 25406, 2681,
                                                                       2747, 11954, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 27980, 0, 3,
                                                                       25406, 10469, 25736, 2747,
                                                                       2813, 12152, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 28376, 0, 3,
                                                                       25736, 10634, 26066, 2813,
                                                                       2879, 12350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 28772, 0, 3,
                                                                       26396, 11360, 26792, 3011,
                                                                       3089, 13016, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 29240, 0, 3,
                                                                       26792, 11558, 27188, 3089,
                                                                       3167, 13250, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 29708, 0, 3,
                                                                       27188, 11756, 27584, 3167,
                                                                       3245, 13484, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 30176, 0, 3,
                                                                       27584, 11954, 27980, 3245,
                                                                       3323, 13718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 30644, 0, 3,
                                                                       27980, 12152, 28376, 3323,
                                                                       3401, 13952, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31112, 3, 3557,
                                                                       3560, 14186, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31122, 3, 3560,
                                                                       3563, 14192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31132, 3, 3563,
                                                                       3566, 14198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31142, 3, 3566,
                                                                       3569, 14204, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31152, 3, 3569,
                                                                       3572, 14210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31162, 3, 3572,
                                                                       3575, 14216, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31172, 3, 3575,
                                                                       3578, 14222, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31182, 3, 3578,
                                                                       3581, 14228, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31192, 3, 3581,
                                                                       3584, 14234, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31202, 3, 3584,
                                                                       3587, 14240, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31212, 3, 3587,
                                                                       3590, 14246, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31222, 3, 3590,
                                                                       3593, 14252, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31232, 3, 3593,
                                                                       3596, 14258, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31242, 3, 3596,
                                                                       3599, 14264, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31252, 3, 3599,
                                                                       3602, 14270, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31262, 3, 3602,
                                                                       3605, 14276, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31272, 0, 3,
                                                                       31112, 14186, 31122,
                                                                       14282, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31302, 0, 3,
                                                                       31122, 14192, 31132,
                                                                       14300, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31332, 0, 3,
                                                                       31132, 14198, 31142,
                                                                       14318, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31362, 0, 3,
                                                                       31142, 14204, 31152,
                                                                       14336, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31392, 0, 3,
                                                                       31152, 14210, 31162,
                                                                       14354, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31422, 0, 3,
                                                                       31162, 14216, 31172,
                                                                       14372, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31452, 0, 3,
                                                                       31172, 14222, 31182,
                                                                       14390, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31482, 0, 3,
                                                                       31182, 14228, 31192,
                                                                       14408, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31512, 0, 3,
                                                                       31192, 14234, 31202,
                                                                       14426, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31542, 0, 3,
                                                                       31202, 14240, 31212,
                                                                       14444, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31572, 0, 3,
                                                                       31212, 14246, 31222,
                                                                       14462, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31602, 0, 3,
                                                                       31222, 14252, 31232,
                                                                       14480, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31632, 0, 3,
                                                                       31232, 14258, 31242,
                                                                       14498, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31662, 0, 3,
                                                                       31242, 14264, 31252,
                                                                       14516, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31692, 0, 3,
                                                                       31252, 14270, 31262,
                                                                       14534, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 31722, 0, 3,
                                                                       31272, 14282, 31302, 3746,
                                                                       3764, 14552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 31782, 0, 3,
                                                                       31302, 14300, 31332, 3764,
                                                                       3782, 14588, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 31842, 0, 3,
                                                                       31332, 14318, 31362, 3782,
                                                                       3800, 14624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 31902, 0, 3,
                                                                       31362, 14336, 31392, 3800,
                                                                       3818, 14660, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 31962, 0, 3,
                                                                       31392, 14354, 31422, 3818,
                                                                       3836, 14696, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32022, 0, 3,
                                                                       31422, 14372, 31452, 3836,
                                                                       3854, 14732, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32082, 0, 3,
                                                                       31452, 14390, 31482, 3854,
                                                                       3872, 14768, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32142, 0, 3,
                                                                       31482, 14408, 31512, 3872,
                                                                       3890, 14804, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32202, 0, 3,
                                                                       31512, 14426, 31542, 3890,
                                                                       3908, 14840, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32262, 0, 3,
                                                                       31542, 14444, 31572, 3908,
                                                                       3926, 14876, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32322, 0, 3,
                                                                       31572, 14462, 31602, 3926,
                                                                       3944, 14912, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32382, 0, 3,
                                                                       31602, 14480, 31632, 3944,
                                                                       3962, 14948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32442, 0, 3,
                                                                       31632, 14498, 31662, 3962,
                                                                       3980, 14984, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32502, 0, 3,
                                                                       31662, 14516, 31692, 3980,
                                                                       3998, 15020, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 32562, 0, 3,
                                                                       31722, 14552, 31782, 4034,
                                                                       4064, 15056, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 32662, 0, 3,
                                                                       31782, 14588, 31842, 4064,
                                                                       4094, 15116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 32762, 0, 3,
                                                                       31842, 14624, 31902, 4094,
                                                                       4124, 15176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 32862, 0, 3,
                                                                       31902, 14660, 31962, 4124,
                                                                       4154, 15236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 32962, 0, 3,
                                                                       31962, 14696, 32022, 4154,
                                                                       4184, 15296, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33062, 0, 3,
                                                                       32022, 14732, 32082, 4184,
                                                                       4214, 15356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33162, 0, 3,
                                                                       32082, 14768, 32142, 4214,
                                                                       4244, 15416, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33262, 0, 3,
                                                                       32142, 14804, 32202, 4244,
                                                                       4274, 15476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33362, 0, 3,
                                                                       32202, 14840, 32262, 4274,
                                                                       4304, 15536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33462, 0, 3,
                                                                       32262, 14876, 32322, 4304,
                                                                       4334, 15596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33562, 0, 3,
                                                                       32322, 14912, 32382, 4334,
                                                                       4364, 15656, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33662, 0, 3,
                                                                       32382, 14948, 32442, 4364,
                                                                       4394, 15716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33762, 0, 3,
                                                                       32442, 14984, 32502, 4394,
                                                                       4424, 15776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 33862, 0, 3,
                                                                       32562, 15056, 32662, 4484,
                                                                       4529, 15836, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34012, 0, 3,
                                                                       32662, 15116, 32762, 4529,
                                                                       4574, 15926, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34162, 0, 3,
                                                                       32762, 15176, 32862, 4574,
                                                                       4619, 16016, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34312, 0, 3,
                                                                       32862, 15236, 32962, 4619,
                                                                       4664, 16106, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34462, 0, 3,
                                                                       32962, 15296, 33062, 4664,
                                                                       4709, 16196, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34612, 0, 3,
                                                                       33062, 15356, 33162, 4709,
                                                                       4754, 16286, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34762, 0, 3,
                                                                       33162, 15416, 33262, 4754,
                                                                       4799, 16376, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34912, 0, 3,
                                                                       33262, 15476, 33362, 4799,
                                                                       4844, 16466, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 35062, 0, 3,
                                                                       33362, 15536, 33462, 4844,
                                                                       4889, 16556, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 35212, 0, 3,
                                                                       33462, 15596, 33562, 4889,
                                                                       4934, 16646, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 35362, 0, 3,
                                                                       33562, 15656, 33662, 4934,
                                                                       4979, 16736, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 35512, 0, 3,
                                                                       33662, 15716, 33762, 4979,
                                                                       5024, 16826, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35662, 0, 3,
                                                                       33862, 15836, 34012, 5114,
                                                                       5177, 16916, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35872, 0, 3,
                                                                       34012, 15926, 34162, 5177,
                                                                       5240, 17042, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 36082, 0, 3,
                                                                       34162, 16016, 34312, 5240,
                                                                       5303, 17168, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 36292, 0, 3,
                                                                       34312, 16106, 34462, 5303,
                                                                       5366, 17294, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 36502, 0, 3,
                                                                       34462, 16196, 34612, 5366,
                                                                       5429, 17420, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 36712, 0, 3,
                                                                       34612, 16286, 34762, 5429,
                                                                       5492, 17546, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 36922, 0, 3,
                                                                       34762, 16376, 34912, 5492,
                                                                       5555, 17672, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 37132, 0, 3,
                                                                       34912, 16466, 35062, 5555,
                                                                       5618, 17798, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 37342, 0, 3,
                                                                       35062, 16556, 35212, 5618,
                                                                       5681, 17924, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 37552, 0, 3,
                                                                       35212, 16646, 35362, 5681,
                                                                       5744, 18050, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 37762, 0, 3,
                                                                       35362, 16736, 35512, 5744,
                                                                       5807, 18176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 37972, 0, 3,
                                                                       35662, 16916, 35872, 5933,
                                                                       6017, 18302, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38252, 0, 3,
                                                                       35872, 17042, 36082, 6017,
                                                                       6101, 18470, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38532, 0, 3,
                                                                       36082, 17168, 36292, 6101,
                                                                       6185, 18638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38812, 0, 3,
                                                                       36292, 17294, 36502, 6185,
                                                                       6269, 18806, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39092, 0, 3,
                                                                       36502, 17420, 36712, 6269,
                                                                       6353, 18974, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39372, 0, 3,
                                                                       36712, 17546, 36922, 6353,
                                                                       6437, 19142, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39652, 0, 3,
                                                                       36922, 17672, 37132, 6437,
                                                                       6521, 19310, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39932, 0, 3,
                                                                       37132, 17798, 37342, 6521,
                                                                       6605, 19478, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 40212, 0, 3,
                                                                       37342, 17924, 37552, 6605,
                                                                       6689, 19646, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 40492, 0, 3,
                                                                       37552, 18050, 37762, 6689,
                                                                       6773, 19814, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 40772, 0, 3,
                                                                       37972, 18302, 38252, 6941,
                                                                       7049, 19982, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41132, 0, 3,
                                                                       38252, 18470, 38532, 7049,
                                                                       7157, 20198, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41492, 0, 3,
                                                                       38532, 18638, 38812, 7157,
                                                                       7265, 20414, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41852, 0, 3,
                                                                       38812, 18806, 39092, 7265,
                                                                       7373, 20630, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42212, 0, 3,
                                                                       39092, 18974, 39372, 7373,
                                                                       7481, 20846, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42572, 0, 3,
                                                                       39372, 19142, 39652, 7481,
                                                                       7589, 21062, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42932, 0, 3,
                                                                       39652, 19310, 39932, 7589,
                                                                       7697, 21278, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 43292, 0, 3,
                                                                       39932, 19478, 40212, 7697,
                                                                       7805, 21494, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 43652, 0, 3,
                                                                       40212, 19646, 40492, 7805,
                                                                       7913, 21710, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 44012, 0, 3,
                                                                       40772, 19982, 41132, 8129,
                                                                       8264, 21926, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 44462, 0, 3,
                                                                       41132, 20198, 41492, 8264,
                                                                       8399, 22196, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 44912, 0, 3,
                                                                       41492, 20414, 41852, 8399,
                                                                       8534, 22466, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 45362, 0, 3,
                                                                       41852, 20630, 42212, 8534,
                                                                       8669, 22736, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 45812, 0, 3,
                                                                       42212, 20846, 42572, 8669,
                                                                       8804, 23006, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 46262, 0, 3,
                                                                       42572, 21062, 42932, 8804,
                                                                       8939, 23276, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 46712, 0, 3,
                                                                       42932, 21278, 43292, 8939,
                                                                       9074, 23546, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 47162, 0, 3,
                                                                       43292, 21494, 43652, 9074,
                                                                       9209, 23816, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 47612, 0, 3,
                                                                       44012, 21926, 44462, 9479,
                                                                       9644, 24086, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 48162, 0, 3,
                                                                       44462, 22196, 44912, 9644,
                                                                       9809, 24416, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 48712, 0, 3,
                                                                       44912, 22466, 45362, 9809,
                                                                       9974, 24746, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 49262, 0, 3,
                                                                       45362, 22736, 45812, 9974,
                                                                       10139, 25076, ncols,
                                                                       gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 49812, 0, 3,
                                                                       45812, 23006, 46262,
                                                                       10139, 10304, 25406,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 50362, 0, 3,
                                                                       46262, 23276, 46712,
                                                                       10304, 10469, 25736,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 50912, 0, 3,
                                                                       46712, 23546, 47162,
                                                                       10469, 10634, 26066,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 51462, 0, 3,
                                                                       47612, 24086, 48162,
                                                                       10964, 11162, 26396,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 52122, 0, 3,
                                                                       48162, 24416, 48712,
                                                                       11162, 11360, 26792,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 52782, 0, 3,
                                                                       48712, 24746, 49262,
                                                                       11360, 11558, 27188,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 53442, 0, 3,
                                                                       49262, 25076, 49812,
                                                                       11558, 11756, 27584,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 54102, 0, 3,
                                                                       49812, 25406, 50362,
                                                                       11756, 11954, 27980,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 54762, 0, 3,
                                                                       50362, 25736, 50912,
                                                                       11954, 12152, 28376,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 55422, 0, 3,
                                                                       51462, 26396, 52122,
                                                                       12548, 12782, 28772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 56202, 0, 3,
                                                                       52122, 26792, 52782,
                                                                       12782, 13016, 29240,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 56982, 0, 3,
                                                                       52782, 27188, 53442,
                                                                       13016, 13250, 29708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 57762, 0, 3,
                                                                       53442, 27584, 54102,
                                                                       13250, 13484, 30176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 58542, 0, 3,
                                                                       54102, 27980, 54762,
                                                                       13484, 13718, 30644,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59322, 3, 14186,
                                                                       14192, 31132, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59337, 3, 14192,
                                                                       14198, 31142, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59352, 3, 14198,
                                                                       14204, 31152, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59367, 3, 14204,
                                                                       14210, 31162, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59382, 3, 14210,
                                                                       14216, 31172, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59397, 3, 14216,
                                                                       14222, 31182, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59412, 3, 14222,
                                                                       14228, 31192, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59427, 3, 14228,
                                                                       14234, 31202, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59442, 3, 14234,
                                                                       14240, 31212, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59457, 3, 14240,
                                                                       14246, 31222, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59472, 3, 14246,
                                                                       14252, 31232, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59487, 3, 14252,
                                                                       14258, 31242, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59502, 3, 14258,
                                                                       14264, 31252, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59517, 3, 14264,
                                                                       14270, 31262, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59532, 0, 3,
                                                                       59322, 31132, 59337,
                                                                       14282, 14300, 31332,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59577, 0, 3,
                                                                       59337, 31142, 59352,
                                                                       14300, 14318, 31362,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59622, 0, 3,
                                                                       59352, 31152, 59367,
                                                                       14318, 14336, 31392,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59667, 0, 3,
                                                                       59367, 31162, 59382,
                                                                       14336, 14354, 31422,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59712, 0, 3,
                                                                       59382, 31172, 59397,
                                                                       14354, 14372, 31452,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59757, 0, 3,
                                                                       59397, 31182, 59412,
                                                                       14372, 14390, 31482,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59802, 0, 3,
                                                                       59412, 31192, 59427,
                                                                       14390, 14408, 31512,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59847, 0, 3,
                                                                       59427, 31202, 59442,
                                                                       14408, 14426, 31542,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59892, 0, 3,
                                                                       59442, 31212, 59457,
                                                                       14426, 14444, 31572,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59937, 0, 3,
                                                                       59457, 31222, 59472,
                                                                       14444, 14462, 31602,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59982, 0, 3,
                                                                       59472, 31232, 59487,
                                                                       14462, 14480, 31632,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 60027, 0, 3,
                                                                       59487, 31242, 59502,
                                                                       14480, 14498, 31662,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 60072, 0, 3,
                                                                       59502, 31252, 59517,
                                                                       14498, 14516, 31692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60117, 0, 3,
                                                                       59532, 31332, 59577,
                                                                       14552, 14588, 31842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60207, 0, 3,
                                                                       59577, 31362, 59622,
                                                                       14588, 14624, 31902,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60297, 0, 3,
                                                                       59622, 31392, 59667,
                                                                       14624, 14660, 31962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60387, 0, 3,
                                                                       59667, 31422, 59712,
                                                                       14660, 14696, 32022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60477, 0, 3,
                                                                       59712, 31452, 59757,
                                                                       14696, 14732, 32082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60567, 0, 3,
                                                                       59757, 31482, 59802,
                                                                       14732, 14768, 32142,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60657, 0, 3,
                                                                       59802, 31512, 59847,
                                                                       14768, 14804, 32202,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60747, 0, 3,
                                                                       59847, 31542, 59892,
                                                                       14804, 14840, 32262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60837, 0, 3,
                                                                       59892, 31572, 59937,
                                                                       14840, 14876, 32322,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60927, 0, 3,
                                                                       59937, 31602, 59982,
                                                                       14876, 14912, 32382,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 61017, 0, 3,
                                                                       59982, 31632, 60027,
                                                                       14912, 14948, 32442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 61107, 0, 3,
                                                                       60027, 31662, 60072,
                                                                       14948, 14984, 32502,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61197, 0, 3,
                                                                       60117, 31842, 60207,
                                                                       15056, 15116, 32762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61347, 0, 3,
                                                                       60207, 31902, 60297,
                                                                       15116, 15176, 32862,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61497, 0, 3,
                                                                       60297, 31962, 60387,
                                                                       15176, 15236, 32962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61647, 0, 3,
                                                                       60387, 32022, 60477,
                                                                       15236, 15296, 33062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61797, 0, 3,
                                                                       60477, 32082, 60567,
                                                                       15296, 15356, 33162,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61947, 0, 3,
                                                                       60567, 32142, 60657,
                                                                       15356, 15416, 33262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 62097, 0, 3,
                                                                       60657, 32202, 60747,
                                                                       15416, 15476, 33362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 62247, 0, 3,
                                                                       60747, 32262, 60837,
                                                                       15476, 15536, 33462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 62397, 0, 3,
                                                                       60837, 32322, 60927,
                                                                       15536, 15596, 33562,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 62547, 0, 3,
                                                                       60927, 32382, 61017,
                                                                       15596, 15656, 33662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 62697, 0, 3,
                                                                       61017, 32442, 61107,
                                                                       15656, 15716, 33762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 62847, 0, 3,
                                                                       61197, 32762, 61347,
                                                                       15836, 15926, 34162,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63072, 0, 3,
                                                                       61347, 32862, 61497,
                                                                       15926, 16016, 34312,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63297, 0, 3,
                                                                       61497, 32962, 61647,
                                                                       16016, 16106, 34462,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63522, 0, 3,
                                                                       61647, 33062, 61797,
                                                                       16106, 16196, 34612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63747, 0, 3,
                                                                       61797, 33162, 61947,
                                                                       16196, 16286, 34762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63972, 0, 3,
                                                                       61947, 33262, 62097,
                                                                       16286, 16376, 34912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64197, 0, 3,
                                                                       62097, 33362, 62247,
                                                                       16376, 16466, 35062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64422, 0, 3,
                                                                       62247, 33462, 62397,
                                                                       16466, 16556, 35212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64647, 0, 3,
                                                                       62397, 33562, 62547,
                                                                       16556, 16646, 35362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64872, 0, 3,
                                                                       62547, 33662, 62697,
                                                                       16646, 16736, 35512,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 65097, 0, 3,
                                                                       62847, 34162, 63072,
                                                                       16916, 17042, 36082,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 65412, 0, 3,
                                                                       63072, 34312, 63297,
                                                                       17042, 17168, 36292,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 65727, 0, 3,
                                                                       63297, 34462, 63522,
                                                                       17168, 17294, 36502,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66042, 0, 3,
                                                                       63522, 34612, 63747,
                                                                       17294, 17420, 36712,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66357, 0, 3,
                                                                       63747, 34762, 63972,
                                                                       17420, 17546, 36922,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66672, 0, 3,
                                                                       63972, 34912, 64197,
                                                                       17546, 17672, 37132,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66987, 0, 3,
                                                                       64197, 35062, 64422,
                                                                       17672, 17798, 37342,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 67302, 0, 3,
                                                                       64422, 35212, 64647,
                                                                       17798, 17924, 37552,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 67617, 0, 3,
                                                                       64647, 35362, 64872,
                                                                       17924, 18050, 37762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 67932, 0, 3,
                                                                       65097, 36082, 65412,
                                                                       18302, 18470, 38532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 68352, 0, 3,
                                                                       65412, 36292, 65727,
                                                                       18470, 18638, 38812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 68772, 0, 3,
                                                                       65727, 36502, 66042,
                                                                       18638, 18806, 39092,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 69192, 0, 3,
                                                                       66042, 36712, 66357,
                                                                       18806, 18974, 39372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 69612, 0, 3,
                                                                       66357, 36922, 66672,
                                                                       18974, 19142, 39652,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 70032, 0, 3,
                                                                       66672, 37132, 66987,
                                                                       19142, 19310, 39932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 70452, 0, 3,
                                                                       66987, 37342, 67302,
                                                                       19310, 19478, 40212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 70872, 0, 3,
                                                                       67302, 37552, 67617,
                                                                       19478, 19646, 40492,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 71292, 0, 3,
                                                                       67932, 38532, 68352,
                                                                       19982, 20198, 41492,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 71832, 0, 3,
                                                                       68352, 38812, 68772,
                                                                       20198, 20414, 41852,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 72372, 0, 3,
                                                                       68772, 39092, 69192,
                                                                       20414, 20630, 42212,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 72912, 0, 3,
                                                                       69192, 39372, 69612,
                                                                       20630, 20846, 42572,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 73452, 0, 3,
                                                                       69612, 39652, 70032,
                                                                       20846, 21062, 42932,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 73992, 0, 3,
                                                                       70032, 39932, 70452,
                                                                       21062, 21278, 43292,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 74532, 0, 3,
                                                                       70452, 40212, 70872,
                                                                       21278, 21494, 43652,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 75072, 0, 3,
                                                                       71292, 41492, 71832,
                                                                       21926, 22196, 44912,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 75747, 0, 3,
                                                                       71832, 41852, 72372,
                                                                       22196, 22466, 45362,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 76422, 0, 3,
                                                                       72372, 42212, 72912,
                                                                       22466, 22736, 45812,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 77097, 0, 3,
                                                                       72912, 42572, 73452,
                                                                       22736, 23006, 46262,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 77772, 0, 3,
                                                                       73452, 42932, 73992,
                                                                       23006, 23276, 46712,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 78447, 0, 3,
                                                                       73992, 43292, 74532,
                                                                       23276, 23546, 47162,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 79122, 0, 3,
                                                                       75072, 44912, 75747,
                                                                       24086, 24416, 48712,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 79947, 0, 3,
                                                                       75747, 45362, 76422,
                                                                       24416, 24746, 49262,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 80772, 0, 3,
                                                                       76422, 45812, 77097,
                                                                       24746, 25076, 49812,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 81597, 0, 3,
                                                                       77097, 46262, 77772,
                                                                       25076, 25406, 50362,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 82422, 0, 3,
                                                                       77772, 46712, 78447,
                                                                       25406, 25736, 50912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 83247, 0, 3,
                                                                       79122, 48712, 79947,
                                                                       26396, 26792, 52782,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 84237, 0, 3,
                                                                       79947, 49262, 80772,
                                                                       26792, 27188, 53442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 85227, 0, 3,
                                                                       80772, 49812, 81597,
                                                                       27188, 27584, 54102,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 86217, 0, 3,
                                                                       81597, 50362, 82422,
                                                                       27584, 27980, 54762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 87207, 0, 3,
                                                                       83247, 52782, 84237,
                                                                       28772, 29240, 56982,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 88377, 0, 3,
                                                                       84237, 53442, 85227,
                                                                       29240, 29708, 57762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 89547, 0, 3,
                                                                       85227, 54102, 86217,
                                                                       29708, 30176, 58542,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90717, 3, 31112,
                                                                       31122, 59322, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90738, 3, 31122,
                                                                       31132, 59337, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90759, 3, 31132,
                                                                       31142, 59352, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90780, 3, 31142,
                                                                       31152, 59367, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90801, 3, 31152,
                                                                       31162, 59382, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90822, 3, 31162,
                                                                       31172, 59397, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90843, 3, 31172,
                                                                       31182, 59412, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90864, 3, 31182,
                                                                       31192, 59427, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90885, 3, 31192,
                                                                       31202, 59442, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90906, 3, 31202,
                                                                       31212, 59457, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90927, 3, 31212,
                                                                       31222, 59472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90948, 3, 31222,
                                                                       31232, 59487, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90969, 3, 31232,
                                                                       31242, 59502, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90990, 3, 31242,
                                                                       31252, 59517, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91011, 0, 3,
                                                                       90717, 59322, 90738,
                                                                       31272, 31302, 59532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91074, 0, 3,
                                                                       90738, 59337, 90759,
                                                                       31302, 31332, 59577,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91137, 0, 3,
                                                                       90759, 59352, 90780,
                                                                       31332, 31362, 59622,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91200, 0, 3,
                                                                       90780, 59367, 90801,
                                                                       31362, 31392, 59667,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91263, 0, 3,
                                                                       90801, 59382, 90822,
                                                                       31392, 31422, 59712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91326, 0, 3,
                                                                       90822, 59397, 90843,
                                                                       31422, 31452, 59757,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91389, 0, 3,
                                                                       90843, 59412, 90864,
                                                                       31452, 31482, 59802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91452, 0, 3,
                                                                       90864, 59427, 90885,
                                                                       31482, 31512, 59847,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91515, 0, 3,
                                                                       90885, 59442, 90906,
                                                                       31512, 31542, 59892,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91578, 0, 3,
                                                                       90906, 59457, 90927,
                                                                       31542, 31572, 59937,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91641, 0, 3,
                                                                       90927, 59472, 90948,
                                                                       31572, 31602, 59982,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91704, 0, 3,
                                                                       90948, 59487, 90969,
                                                                       31602, 31632, 60027,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91767, 0, 3,
                                                                       90969, 59502, 90990,
                                                                       31632, 31662, 60072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 91830, 0, 3,
                                                                       91011, 59532, 91074,
                                                                       31722, 31782, 60117,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 91956, 0, 3,
                                                                       91074, 59577, 91137,
                                                                       31782, 31842, 60207,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92082, 0, 3,
                                                                       91137, 59622, 91200,
                                                                       31842, 31902, 60297,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92208, 0, 3,
                                                                       91200, 59667, 91263,
                                                                       31902, 31962, 60387,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92334, 0, 3,
                                                                       91263, 59712, 91326,
                                                                       31962, 32022, 60477,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92460, 0, 3,
                                                                       91326, 59757, 91389,
                                                                       32022, 32082, 60567,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92586, 0, 3,
                                                                       91389, 59802, 91452,
                                                                       32082, 32142, 60657,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92712, 0, 3,
                                                                       91452, 59847, 91515,
                                                                       32142, 32202, 60747,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92838, 0, 3,
                                                                       91515, 59892, 91578,
                                                                       32202, 32262, 60837,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92964, 0, 3,
                                                                       91578, 59937, 91641,
                                                                       32262, 32322, 60927,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93090, 0, 3,
                                                                       91641, 59982, 91704,
                                                                       32322, 32382, 61017,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93216, 0, 3,
                                                                       91704, 60027, 91767,
                                                                       32382, 32442, 61107,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 93342, 0, 3,
                                                                       91830, 60117, 91956,
                                                                       32562, 32662, 61197,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 93552, 0, 3,
                                                                       91956, 60207, 92082,
                                                                       32662, 32762, 61347,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 93762, 0, 3,
                                                                       92082, 60297, 92208,
                                                                       32762, 32862, 61497,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 93972, 0, 3,
                                                                       92208, 60387, 92334,
                                                                       32862, 32962, 61647,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94182, 0, 3,
                                                                       92334, 60477, 92460,
                                                                       32962, 33062, 61797,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94392, 0, 3,
                                                                       92460, 60567, 92586,
                                                                       33062, 33162, 61947,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94602, 0, 3,
                                                                       92586, 60657, 92712,
                                                                       33162, 33262, 62097,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94812, 0, 3,
                                                                       92712, 60747, 92838,
                                                                       33262, 33362, 62247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95022, 0, 3,
                                                                       92838, 60837, 92964,
                                                                       33362, 33462, 62397,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95232, 0, 3,
                                                                       92964, 60927, 93090,
                                                                       33462, 33562, 62547,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95442, 0, 3,
                                                                       93090, 61017, 93216,
                                                                       33562, 33662, 62697,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 95652, 0, 3,
                                                                       93342, 61197, 93552,
                                                                       33862, 34012, 62847,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 95967, 0, 3,
                                                                       93552, 61347, 93762,
                                                                       34012, 34162, 63072,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 96282, 0, 3,
                                                                       93762, 61497, 93972,
                                                                       34162, 34312, 63297,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 96597, 0, 3,
                                                                       93972, 61647, 94182,
                                                                       34312, 34462, 63522,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 96912, 0, 3,
                                                                       94182, 61797, 94392,
                                                                       34462, 34612, 63747,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 97227, 0, 3,
                                                                       94392, 61947, 94602,
                                                                       34612, 34762, 63972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 97542, 0, 3,
                                                                       94602, 62097, 94812,
                                                                       34762, 34912, 64197,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 97857, 0, 3,
                                                                       94812, 62247, 95022,
                                                                       34912, 35062, 64422,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 98172, 0, 3,
                                                                       95022, 62397, 95232,
                                                                       35062, 35212, 64647,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 98487, 0, 3,
                                                                       95232, 62547, 95442,
                                                                       35212, 35362, 64872,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 98802, 0, 3,
                                                                       95652, 62847, 95967,
                                                                       35662, 35872, 65097,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 99243, 0, 3,
                                                                       95967, 63072, 96282,
                                                                       35872, 36082, 65412,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 99684, 0, 3,
                                                                       96282, 63297, 96597,
                                                                       36082, 36292, 65727,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 100125, 0, 3,
                                                                       96597, 63522, 96912,
                                                                       36292, 36502, 66042,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 100566, 0, 3,
                                                                       96912, 63747, 97227,
                                                                       36502, 36712, 66357,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 101007, 0, 3,
                                                                       97227, 63972, 97542,
                                                                       36712, 36922, 66672,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 101448, 0, 3,
                                                                       97542, 64197, 97857,
                                                                       36922, 37132, 66987,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 101889, 0, 3,
                                                                       97857, 64422, 98172,
                                                                       37132, 37342, 67302,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 102330, 0, 3,
                                                                       98172, 64647, 98487,
                                                                       37342, 37552, 67617,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 102771, 0, 3,
                                                                       98802, 65097, 99243,
                                                                       37972, 38252, 67932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 103359, 0, 3,
                                                                       99243, 65412, 99684,
                                                                       38252, 38532, 68352,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 103947, 0, 3,
                                                                       99684, 65727, 100125,
                                                                       38532, 38812, 68772,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 104535, 0, 3,
                                                                       100125, 66042, 100566,
                                                                       38812, 39092, 69192,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 105123, 0, 3,
                                                                       100566, 66357, 101007,
                                                                       39092, 39372, 69612,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 105711, 0, 3,
                                                                       101007, 66672, 101448,
                                                                       39372, 39652, 70032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 106299, 0, 3,
                                                                       101448, 66987, 101889,
                                                                       39652, 39932, 70452,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 106887, 0, 3,
                                                                       101889, 67302, 102330,
                                                                       39932, 40212, 70872,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 107475, 0, 3,
                                                                       102771, 67932, 103359,
                                                                       40772, 41132, 71292,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 108231, 0, 3,
                                                                       103359, 68352, 103947,
                                                                       41132, 41492, 71832,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 108987, 0, 3,
                                                                       103947, 68772, 104535,
                                                                       41492, 41852, 72372,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 109743, 0, 3,
                                                                       104535, 69192, 105123,
                                                                       41852, 42212, 72912,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 110499, 0, 3,
                                                                       105123, 69612, 105711,
                                                                       42212, 42572, 73452,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 111255, 0, 3,
                                                                       105711, 70032, 106299,
                                                                       42572, 42932, 73992,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 112011, 0, 3,
                                                                       106299, 70452, 106887,
                                                                       42932, 43292, 74532,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 112767, 0, 3,
                                                                       107475, 71292, 108231,
                                                                       44012, 44462, 75072,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 113712, 0, 3,
                                                                       108231, 71832, 108987,
                                                                       44462, 44912, 75747,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 114657, 0, 3,
                                                                       108987, 72372, 109743,
                                                                       44912, 45362, 76422,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 115602, 0, 3,
                                                                       109743, 72912, 110499,
                                                                       45362, 45812, 77097,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 116547, 0, 3,
                                                                       110499, 73452, 111255,
                                                                       45812, 46262, 77772,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 117492, 0, 3,
                                                                       111255, 73992, 112011,
                                                                       46262, 46712, 78447,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 118437, 0, 3,
                                                                       112767, 75072, 113712,
                                                                       47612, 48162, 79122,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 119592, 0, 3,
                                                                       113712, 75747, 114657,
                                                                       48162, 48712, 79947,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 120747, 0, 3,
                                                                       114657, 76422, 115602,
                                                                       48712, 49262, 80772,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 121902, 0, 3,
                                                                       115602, 77097, 116547,
                                                                       49262, 49812, 81597,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 123057, 0, 3,
                                                                       116547, 77772, 117492,
                                                                       49812, 50362, 82422,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 124212, 0, 3,
                                                                       118437, 79122, 119592,
                                                                       51462, 52122, 83247,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 125598, 0, 3,
                                                                       119592, 79947, 120747,
                                                                       52122, 52782, 84237,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 126984, 0, 3,
                                                                       120747, 80772, 121902,
                                                                       52782, 53442, 85227,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 128370, 0, 3,
                                                                       121902, 81597, 123057,
                                                                       53442, 54102, 86217,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 129756, 0, 3,
                                                                       124212, 83247, 125598,
                                                                       55422, 56202, 87207,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 131394, 0, 3,
                                                                       125598, 84237, 126984,
                                                                       56202, 56982, 88377,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 133032, 0, 3,
                                                                       126984, 85227, 128370,
                                                                       56982, 57762, 89547,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134670, 3, 59322,
                                                                       59337, 90759, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134698, 3, 59337,
                                                                       59352, 90780, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134726, 3, 59352,
                                                                       59367, 90801, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134754, 3, 59367,
                                                                       59382, 90822, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134782, 3, 59382,
                                                                       59397, 90843, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134810, 3, 59397,
                                                                       59412, 90864, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134838, 3, 59412,
                                                                       59427, 90885, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134866, 3, 59427,
                                                                       59442, 90906, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134894, 3, 59442,
                                                                       59457, 90927, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134922, 3, 59457,
                                                                       59472, 90948, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134950, 3, 59472,
                                                                       59487, 90969, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134978, 3, 59487,
                                                                       59502, 90990, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135006, 0, 3,
                                                                       134670, 90759, 134698,
                                                                       59532, 59577, 91137,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135090, 0, 3,
                                                                       134698, 90780, 134726,
                                                                       59577, 59622, 91200,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135174, 0, 3,
                                                                       134726, 90801, 134754,
                                                                       59622, 59667, 91263,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135258, 0, 3,
                                                                       134754, 90822, 134782,
                                                                       59667, 59712, 91326,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135342, 0, 3,
                                                                       134782, 90843, 134810,
                                                                       59712, 59757, 91389,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135426, 0, 3,
                                                                       134810, 90864, 134838,
                                                                       59757, 59802, 91452,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135510, 0, 3,
                                                                       134838, 90885, 134866,
                                                                       59802, 59847, 91515,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135594, 0, 3,
                                                                       134866, 90906, 134894,
                                                                       59847, 59892, 91578,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135678, 0, 3,
                                                                       134894, 90927, 134922,
                                                                       59892, 59937, 91641,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135762, 0, 3,
                                                                       134922, 90948, 134950,
                                                                       59937, 59982, 91704,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135846, 0, 3,
                                                                       134950, 90969, 134978,
                                                                       59982, 60027, 91767,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 135930, 0, 3,
                                                                       135006, 91137, 135090,
                                                                       60117, 60207, 92082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136098, 0, 3,
                                                                       135090, 91200, 135174,
                                                                       60207, 60297, 92208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136266, 0, 3,
                                                                       135174, 91263, 135258,
                                                                       60297, 60387, 92334,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136434, 0, 3,
                                                                       135258, 91326, 135342,
                                                                       60387, 60477, 92460,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136602, 0, 3,
                                                                       135342, 91389, 135426,
                                                                       60477, 60567, 92586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136770, 0, 3,
                                                                       135426, 91452, 135510,
                                                                       60567, 60657, 92712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136938, 0, 3,
                                                                       135510, 91515, 135594,
                                                                       60657, 60747, 92838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 137106, 0, 3,
                                                                       135594, 91578, 135678,
                                                                       60747, 60837, 92964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 137274, 0, 3,
                                                                       135678, 91641, 135762,
                                                                       60837, 60927, 93090,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 137442, 0, 3,
                                                                       135762, 91704, 135846,
                                                                       60927, 61017, 93216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 137610, 0, 3,
                                                                       135930, 92082, 136098,
                                                                       61197, 61347, 93762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 137890, 0, 3,
                                                                       136098, 92208, 136266,
                                                                       61347, 61497, 93972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 138170, 0, 3,
                                                                       136266, 92334, 136434,
                                                                       61497, 61647, 94182,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 138450, 0, 3,
                                                                       136434, 92460, 136602,
                                                                       61647, 61797, 94392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 138730, 0, 3,
                                                                       136602, 92586, 136770,
                                                                       61797, 61947, 94602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 139010, 0, 3,
                                                                       136770, 92712, 136938,
                                                                       61947, 62097, 94812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 139290, 0, 3,
                                                                       136938, 92838, 137106,
                                                                       62097, 62247, 95022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 139570, 0, 3,
                                                                       137106, 92964, 137274,
                                                                       62247, 62397, 95232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 139850, 0, 3,
                                                                       137274, 93090, 137442,
                                                                       62397, 62547, 95442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 140130, 0, 3,
                                                                       137610, 93762, 137890,
                                                                       62847, 63072, 96282,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 140550, 0, 3,
                                                                       137890, 93972, 138170,
                                                                       63072, 63297, 96597,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 140970, 0, 3,
                                                                       138170, 94182, 138450,
                                                                       63297, 63522, 96912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 141390, 0, 3,
                                                                       138450, 94392, 138730,
                                                                       63522, 63747, 97227,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 141810, 0, 3,
                                                                       138730, 94602, 139010,
                                                                       63747, 63972, 97542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 142230, 0, 3,
                                                                       139010, 94812, 139290,
                                                                       63972, 64197, 97857,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 142650, 0, 3,
                                                                       139290, 95022, 139570,
                                                                       64197, 64422, 98172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 143070, 0, 3,
                                                                       139570, 95232, 139850,
                                                                       64422, 64647, 98487,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 143490, 0, 3,
                                                                       140130, 96282, 140550,
                                                                       65097, 65412, 99684,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 144078, 0, 3,
                                                                       140550, 96597, 140970,
                                                                       65412, 65727, 100125,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 144666, 0, 3,
                                                                       140970, 96912, 141390,
                                                                       65727, 66042, 100566,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 145254, 0, 3,
                                                                       141390, 97227, 141810,
                                                                       66042, 66357, 101007,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 145842, 0, 3,
                                                                       141810, 97542, 142230,
                                                                       66357, 66672, 101448,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 146430, 0, 3,
                                                                       142230, 97857, 142650,
                                                                       66672, 66987, 101889,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 147018, 0, 3,
                                                                       142650, 98172, 143070,
                                                                       66987, 67302, 102330,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 147606, 0, 3,
                                                                       143490, 99684, 144078,
                                                                       67932, 68352, 103947,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 148390, 0, 3,
                                                                       144078, 100125, 144666,
                                                                       68352, 68772, 104535,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 149174, 0, 3,
                                                                       144666, 100566, 145254,
                                                                       68772, 69192, 105123,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 149958, 0, 3,
                                                                       145254, 101007, 145842,
                                                                       69192, 69612, 105711,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 150742, 0, 3,
                                                                       145842, 101448, 146430,
                                                                       69612, 70032, 106299,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 151526, 0, 3,
                                                                       146430, 101889, 147018,
                                                                       70032, 70452, 106887,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 152310, 0, 3,
                                                                       147606, 103947, 148390,
                                                                       71292, 71832, 108987,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 153318, 0, 3,
                                                                       148390, 104535, 149174,
                                                                       71832, 72372, 109743,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 154326, 0, 3,
                                                                       149174, 105123, 149958,
                                                                       72372, 72912, 110499,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 155334, 0, 3,
                                                                       149958, 105711, 150742,
                                                                       72912, 73452, 111255,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 156342, 0, 3,
                                                                       150742, 106299, 151526,
                                                                       73452, 73992, 112011,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 157350, 0, 3,
                                                                       152310, 108987, 153318,
                                                                       75072, 75747, 114657,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 158610, 0, 3,
                                                                       153318, 109743, 154326,
                                                                       75747, 76422, 115602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 159870, 0, 3,
                                                                       154326, 110499, 155334,
                                                                       76422, 77097, 116547,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 161130, 0, 3,
                                                                       155334, 111255, 156342,
                                                                       77097, 77772, 117492,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 162390, 0, 3,
                                                                       157350, 114657, 158610,
                                                                       79122, 79947, 120747,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 163930, 0, 3,
                                                                       158610, 115602, 159870,
                                                                       79947, 80772, 121902,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 165470, 0, 3,
                                                                       159870, 116547, 161130,
                                                                       80772, 81597, 123057,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 167010, 0, 3,
                                                                       162390, 120747, 163930,
                                                                       83247, 84237, 126984,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 168858, 0, 3,
                                                                       163930, 121902, 165470,
                                                                       84237, 85227, 128370,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 170706, 0, 3,
                                                                       167010, 126984, 168858,
                                                                       87207, 88377, 133032,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172890, 3, 90717,
                                                                       90738, 134670, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172926, 3, 90738,
                                                                       90759, 134698, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172962, 3, 90759,
                                                                       90780, 134726, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172998, 3, 90780,
                                                                       90801, 134754, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173034, 3, 90801,
                                                                       90822, 134782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173070, 3, 90822,
                                                                       90843, 134810, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173106, 3, 90843,
                                                                       90864, 134838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173142, 3, 90864,
                                                                       90885, 134866, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173178, 3, 90885,
                                                                       90906, 134894, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173214, 3, 90906,
                                                                       90927, 134922, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173250, 3, 90927,
                                                                       90948, 134950, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173286, 3, 90948,
                                                                       90969, 134978, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173322, 0, 3,
                                                                       172890, 134670, 172926,
                                                                       91011, 91074, 135006,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173430, 0, 3,
                                                                       172926, 134698, 172962,
                                                                       91074, 91137, 135090,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173538, 0, 3,
                                                                       172962, 134726, 172998,
                                                                       91137, 91200, 135174,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173646, 0, 3,
                                                                       172998, 134754, 173034,
                                                                       91200, 91263, 135258,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173754, 0, 3,
                                                                       173034, 134782, 173070,
                                                                       91263, 91326, 135342,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173862, 0, 3,
                                                                       173070, 134810, 173106,
                                                                       91326, 91389, 135426,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173970, 0, 3,
                                                                       173106, 134838, 173142,
                                                                       91389, 91452, 135510,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 174078, 0, 3,
                                                                       173142, 134866, 173178,
                                                                       91452, 91515, 135594,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 174186, 0, 3,
                                                                       173178, 134894, 173214,
                                                                       91515, 91578, 135678,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 174294, 0, 3,
                                                                       173214, 134922, 173250,
                                                                       91578, 91641, 135762,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 174402, 0, 3,
                                                                       173250, 134950, 173286,
                                                                       91641, 91704, 135846,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 174510, 0, 3,
                                                                       173322, 135006, 173430,
                                                                       91830, 91956, 135930,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 174726, 0, 3,
                                                                       173430, 135090, 173538,
                                                                       91956, 92082, 136098,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 174942, 0, 3,
                                                                       173538, 135174, 173646,
                                                                       92082, 92208, 136266,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 175158, 0, 3,
                                                                       173646, 135258, 173754,
                                                                       92208, 92334, 136434,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 175374, 0, 3,
                                                                       173754, 135342, 173862,
                                                                       92334, 92460, 136602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 175590, 0, 3,
                                                                       173862, 135426, 173970,
                                                                       92460, 92586, 136770,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 175806, 0, 3,
                                                                       173970, 135510, 174078,
                                                                       92586, 92712, 136938,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 176022, 0, 3,
                                                                       174078, 135594, 174186,
                                                                       92712, 92838, 137106,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 176238, 0, 3,
                                                                       174186, 135678, 174294,
                                                                       92838, 92964, 137274,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 176454, 0, 3,
                                                                       174294, 135762, 174402,
                                                                       92964, 93090, 137442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 176670, 0, 3,
                                                                       174510, 135930, 174726,
                                                                       93342, 93552, 137610,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 177030, 0, 3,
                                                                       174726, 136098, 174942,
                                                                       93552, 93762, 137890,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 177390, 0, 3,
                                                                       174942, 136266, 175158,
                                                                       93762, 93972, 138170,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 177750, 0, 3,
                                                                       175158, 136434, 175374,
                                                                       93972, 94182, 138450,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 178110, 0, 3,
                                                                       175374, 136602, 175590,
                                                                       94182, 94392, 138730,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 178470, 0, 3,
                                                                       175590, 136770, 175806,
                                                                       94392, 94602, 139010,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 178830, 0, 3,
                                                                       175806, 136938, 176022,
                                                                       94602, 94812, 139290,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 179190, 0, 3,
                                                                       176022, 137106, 176238,
                                                                       94812, 95022, 139570,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 179550, 0, 3,
                                                                       176238, 137274, 176454,
                                                                       95022, 95232, 139850,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 179910, 0, 3,
                                                                       176670, 137610, 177030,
                                                                       95652, 95967, 140130,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 180450, 0, 3,
                                                                       177030, 137890, 177390,
                                                                       95967, 96282, 140550,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 180990, 0, 3,
                                                                       177390, 138170, 177750,
                                                                       96282, 96597, 140970,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 181530, 0, 3,
                                                                       177750, 138450, 178110,
                                                                       96597, 96912, 141390,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 182070, 0, 3,
                                                                       178110, 138730, 178470,
                                                                       96912, 97227, 141810,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 182610, 0, 3,
                                                                       178470, 139010, 178830,
                                                                       97227, 97542, 142230,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 183150, 0, 3,
                                                                       178830, 139290, 179190,
                                                                       97542, 97857, 142650,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 183690, 0, 3,
                                                                       179190, 139570, 179550,
                                                                       97857, 98172, 143070,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 184230, 0, 3,
                                                                       179910, 140130, 180450,
                                                                       98802, 99243, 143490,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 184986, 0, 3,
                                                                       180450, 140550, 180990,
                                                                       99243, 99684, 144078,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 185742, 0, 3,
                                                                       180990, 140970, 181530,
                                                                       99684, 100125, 144666,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 186498, 0, 3,
                                                                       181530, 141390, 182070,
                                                                       100125, 100566, 145254,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 187254, 0, 3,
                                                                       182070, 141810, 182610,
                                                                       100566, 101007, 145842,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 188010, 0, 3,
                                                                       182610, 142230, 183150,
                                                                       101007, 101448, 146430,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 188766, 0, 3,
                                                                       183150, 142650, 183690,
                                                                       101448, 101889, 147018,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 189522, 0, 3,
                                                                       184230, 143490, 184986,
                                                                       102771, 103359, 147606,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 190530, 0, 3,
                                                                       184986, 144078, 185742,
                                                                       103359, 103947, 148390,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 191538, 0, 3,
                                                                       185742, 144666, 186498,
                                                                       103947, 104535, 149174,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 192546, 0, 3,
                                                                       186498, 145254, 187254,
                                                                       104535, 105123, 149958,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 193554, 0, 3,
                                                                       187254, 145842, 188010,
                                                                       105123, 105711, 150742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 194562, 0, 3,
                                                                       188010, 146430, 188766,
                                                                       105711, 106299, 151526,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 195570, 0, 3,
                                                                       189522, 147606, 190530,
                                                                       107475, 108231, 152310,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 196866, 0, 3,
                                                                       190530, 148390, 191538,
                                                                       108231, 108987, 153318,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 198162, 0, 3,
                                                                       191538, 149174, 192546,
                                                                       108987, 109743, 154326,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 199458, 0, 3,
                                                                       192546, 149958, 193554,
                                                                       109743, 110499, 155334,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 200754, 0, 3,
                                                                       193554, 150742, 194562,
                                                                       110499, 111255, 156342,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 202050, 0, 3,
                                                                       195570, 152310, 196866,
                                                                       112767, 113712, 157350,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 203670, 0, 3,
                                                                       196866, 153318, 198162,
                                                                       113712, 114657, 158610,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 205290, 0, 3,
                                                                       198162, 154326, 199458,
                                                                       114657, 115602, 159870,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 206910, 0, 3,
                                                                       199458, 155334, 200754,
                                                                       115602, 116547, 161130,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 208530, 0, 3,
                                                                       202050, 157350, 203670,
                                                                       118437, 119592, 162390,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 210510, 0, 3,
                                                                       203670, 158610, 205290,
                                                                       119592, 120747, 163930,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 212490, 0, 3,
                                                                       205290, 159870, 206910,
                                                                       120747, 121902, 165470,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 214470, 0, 3,
                                                                       208530, 162390, 210510,
                                                                       124212, 125598, 167010,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 216846, 0, 3,
                                                                       210510, 163930, 212490,
                                                                       125598, 126984, 168858,
                                                                       ncols, gamma, p, q);

                    compute_prim_sok_three_center_electron_repulsion_0(buffer, 219222, 0, 3,
                                                                       214470, 167010, 216846,
                                                                       129756, 131394, 170706,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 222030, 189522, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 223458, 195570, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 225294, 202050, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 227589, 208530, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 230394, 214470, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 233760, 219222, 2808, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 223038, 222030, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 224754, 223458, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 226914, 225294, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 229569, 227589, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 232770, 230394, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 236568, 233760, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 237738, 223038, 224754, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 238998, 224754, 226914, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 240618, 226914, 229569, 15, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 242643, 229569, 232770, 15, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 245118, 232770, 236568, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 248088, 237738, 238998, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 250608, 238998, 240618, 15, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 253848, 240618, 242643, 15, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 257898, 242643, 245118, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 262848, 248088, 250608, 15, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 267048, 250608, 253848, 15, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 272448, 253848, 257898, 15, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 279198, 262848, 267048, 15, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 285498, 267048, 272448, 15, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 293598, 279198, 285498, 15, nmax);

        simdtrf::transform_i_inner(buffer, 302418, 293598, 21, 15, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 302418, 195, nmax);
    }

    for (size_t m = 0; m < 2145; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
