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

    const auto nmax = simdfunc::prepare_buffer(buffer, 306512, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 306512, 222029, 14538, dimensions);

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

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 6, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                                                        16, 17, 18}, ncols, fj, mu, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 76, 0, 3, 7, 8,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 82, 0, 3, 8, 9,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 88, 0, 3, 9, 10,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 94, 0, 3, 10, 11,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 100, 0, 3, 11, 12,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 106, 0, 3, 12, 13,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 112, 0, 3, 13, 14,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 118, 0, 3, 14, 15,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 124, 0, 3, 15, 16,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 130, 0, 3, 16, 17,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 136, 0, 3, 17, 18,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 142, 0, 3, 18, 19,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 148, 0, 3, 19, 20,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 154, 0, 3, 20, 21,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 160, 0, 3, 21, 22,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 166, 0, 3, 22, 23,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 172, 0, 3, 25, 28,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 28, 31,
                                                                       82, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 31, 34,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 202, 0, 3, 34, 37,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 212, 0, 3, 37, 40,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 40, 43,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 232, 0, 3, 43, 46,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 242, 0, 3, 46, 49,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 252, 0, 3, 49, 52,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 262, 0, 3, 52, 55,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 272, 0, 3, 55, 58,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 282, 0, 3, 58, 61,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 292, 0, 3, 61, 64,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 302, 0, 3, 64, 67,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 312, 0, 3, 67, 70,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 322, 0, 3, 76, 82,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 337, 0, 3, 82, 88,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 352, 0, 3, 88, 94,
                                                                       192, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 367, 0, 3, 94,
                                                                       100, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 382, 0, 3, 100,
                                                                       106, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 397, 0, 3, 106,
                                                                       112, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 412, 0, 3, 112,
                                                                       118, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 427, 0, 3, 118,
                                                                       124, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 442, 0, 3, 124,
                                                                       130, 252, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 457, 0, 3, 130,
                                                                       136, 262, 272, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 472, 0, 3, 136,
                                                                       142, 272, 282, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 487, 0, 3, 142,
                                                                       148, 282, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 502, 0, 3, 148,
                                                                       154, 292, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 517, 0, 3, 154,
                                                                       160, 302, 312, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 532, 0, 3, 172,
                                                                       182, 322, 337, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 553, 0, 3, 182,
                                                                       192, 337, 352, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 574, 0, 3, 192,
                                                                       202, 352, 367, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 595, 0, 3, 202,
                                                                       212, 367, 382, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 616, 0, 3, 212,
                                                                       222, 382, 397, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 637, 0, 3, 222,
                                                                       232, 397, 412, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 658, 0, 3, 232,
                                                                       242, 412, 427, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 679, 0, 3, 242,
                                                                       252, 427, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 700, 0, 3, 252,
                                                                       262, 442, 457, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 721, 0, 3, 262,
                                                                       272, 457, 472, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 742, 0, 3, 272,
                                                                       282, 472, 487, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 763, 0, 3, 282,
                                                                       292, 487, 502, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 784, 0, 3, 292,
                                                                       302, 502, 517, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 805, 0, 3, 322,
                                                                       337, 532, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 833, 0, 3, 337,
                                                                       352, 553, 574, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 861, 0, 3, 352,
                                                                       367, 574, 595, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 889, 0, 3, 367,
                                                                       382, 595, 616, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 917, 0, 3, 382,
                                                                       397, 616, 637, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 945, 0, 3, 397,
                                                                       412, 637, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 973, 0, 3, 412,
                                                                       427, 658, 679, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1001, 0, 3, 427,
                                                                       442, 679, 700, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1029, 0, 3, 442,
                                                                       457, 700, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1057, 0, 3, 457,
                                                                       472, 721, 742, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1085, 0, 3, 472,
                                                                       487, 742, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1113, 0, 3, 487,
                                                                       502, 763, 784, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1141, 0, 3, 532,
                                                                       553, 805, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1177, 0, 3, 553,
                                                                       574, 833, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1213, 0, 3, 574,
                                                                       595, 861, 889, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1249, 0, 3, 595,
                                                                       616, 889, 917, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1285, 0, 3, 616,
                                                                       637, 917, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1321, 0, 3, 637,
                                                                       658, 945, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1357, 0, 3, 658,
                                                                       679, 973, 1001, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1393, 0, 3, 679,
                                                                       700, 1001, 1029, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1429, 0, 3, 700,
                                                                       721, 1029, 1057, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1465, 0, 3, 721,
                                                                       742, 1057, 1085, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1501, 0, 3, 742,
                                                                       763, 1085, 1113, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1537, 0, 3, 805,
                                                                       833, 1141, 1177, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1582, 0, 3, 833,
                                                                       861, 1177, 1213, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1627, 0, 3, 861,
                                                                       889, 1213, 1249, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1672, 0, 3, 889,
                                                                       917, 1249, 1285, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1717, 0, 3, 917,
                                                                       945, 1285, 1321, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1762, 0, 3, 945,
                                                                       973, 1321, 1357, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1807, 0, 3, 973,
                                                                       1001, 1357, 1393, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1852, 0, 3, 1001,
                                                                       1029, 1393, 1429, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1897, 0, 3, 1029,
                                                                       1057, 1429, 1465, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1942, 0, 3, 1057,
                                                                       1085, 1465, 1501, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1987, 0, 3, 1141,
                                                                       1177, 1537, 1582, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2042, 0, 3, 1177,
                                                                       1213, 1582, 1627, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2097, 0, 3, 1213,
                                                                       1249, 1627, 1672, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2152, 0, 3, 1249,
                                                                       1285, 1672, 1717, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2207, 0, 3, 1285,
                                                                       1321, 1717, 1762, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2262, 0, 3, 1321,
                                                                       1357, 1762, 1807, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2317, 0, 3, 1357,
                                                                       1393, 1807, 1852, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2372, 0, 3, 1393,
                                                                       1429, 1852, 1897, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2427, 0, 3, 1429,
                                                                       1465, 1897, 1942, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2482, 0, 3, 1537,
                                                                       1582, 1987, 2042, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2548, 0, 3, 1582,
                                                                       1627, 2042, 2097, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2614, 0, 3, 1627,
                                                                       1672, 2097, 2152, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2680, 0, 3, 1672,
                                                                       1717, 2152, 2207, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2746, 0, 3, 1717,
                                                                       1762, 2207, 2262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2812, 0, 3, 1762,
                                                                       1807, 2262, 2317, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2878, 0, 3, 1807,
                                                                       1852, 2317, 2372, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2944, 0, 3, 1852,
                                                                       1897, 2372, 2427, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3010, 0, 3, 1987,
                                                                       2042, 2482, 2548, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 2042,
                                                                       2097, 2548, 2614, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3166, 0, 3, 2097,
                                                                       2152, 2614, 2680, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3244, 0, 3, 2152,
                                                                       2207, 2680, 2746, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3322, 0, 3, 2207,
                                                                       2262, 2746, 2812, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3400, 0, 3, 2262,
                                                                       2317, 2812, 2878, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3478, 0, 3, 2317,
                                                                       2372, 2878, 2944, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3556, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3559, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3562, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3565, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3568, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3571, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3574, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3577, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3580, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3583, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3586, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3589, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3592, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3595, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3598, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3601, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3604, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3607, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3610, 3, 9, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3619, 3, 10, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3628, 3, 11, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3637, 3, 12, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3646, 3, 13, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3655, 3, 14, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3664, 3, 15, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3673, 3, 16, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3682, 3, 17, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3691, 3, 18, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3700, 3, 19, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3709, 3, 20, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3718, 3, 21, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3727, 3, 22, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3736, 3, 23, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3745, 3, 25, 76,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3763, 3, 28, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3781, 3, 31, 88,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3799, 3, 34, 94,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3817, 3, 37, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3835, 3, 40, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3853, 3, 43, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3871, 3, 46, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3889, 3, 49, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3907, 3, 52, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3925, 3, 55, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3943, 3, 58, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3961, 3, 61, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3979, 3, 64, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3997, 3, 67, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4015, 3, 70, 166,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4033, 3, 76, 172,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4063, 3, 82, 182,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4093, 3, 88, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4123, 3, 94, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4153, 3, 100, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4183, 3, 106, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4213, 3, 112, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4243, 3, 118, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4273, 3, 124, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4303, 3, 130, 262,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4333, 3, 136, 272,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4363, 3, 142, 282,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4393, 3, 148, 292,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4423, 3, 154, 302,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4453, 3, 160, 312,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4483, 3, 172, 322,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4528, 3, 182, 337,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4573, 3, 192, 352,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4618, 3, 202, 367,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4663, 3, 212, 382,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4708, 3, 222, 397,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4753, 3, 232, 412,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4798, 3, 242, 427,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4843, 3, 252, 442,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4888, 3, 262, 457,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4933, 3, 272, 472,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4978, 3, 282, 487,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5023, 3, 292, 502,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5068, 3, 302, 517,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5113, 3, 322, 532,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5176, 3, 337, 553,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5239, 3, 352, 574,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5302, 3, 367, 595,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5365, 3, 382, 616,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5428, 3, 397, 637,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5491, 3, 412, 658,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5554, 3, 427, 679,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5617, 3, 442, 700,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5680, 3, 457, 721,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5743, 3, 472, 742,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5806, 3, 487, 763,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5869, 3, 502, 784,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5932, 3, 532, 805,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6016, 3, 553, 833,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6100, 3, 574, 861,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6184, 3, 595, 889,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6268, 3, 616, 917,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6352, 3, 637, 945,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6436, 3, 658, 973,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6520, 3, 679,
                                                                       1001, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6604, 3, 700,
                                                                       1029, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6688, 3, 721,
                                                                       1057, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6772, 3, 742,
                                                                       1085, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6856, 3, 763,
                                                                       1113, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6940, 3, 805,
                                                                       1141, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7048, 3, 833,
                                                                       1177, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7156, 3, 861,
                                                                       1213, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7264, 3, 889,
                                                                       1249, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7372, 3, 917,
                                                                       1285, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7480, 3, 945,
                                                                       1321, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7588, 3, 973,
                                                                       1357, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7696, 3, 1001,
                                                                       1393, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7804, 3, 1029,
                                                                       1429, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7912, 3, 1057,
                                                                       1465, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8020, 3, 1085,
                                                                       1501, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8128, 3, 1141,
                                                                       1537, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8263, 3, 1177,
                                                                       1582, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8398, 3, 1213,
                                                                       1627, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8533, 3, 1249,
                                                                       1672, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8668, 3, 1285,
                                                                       1717, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8803, 3, 1321,
                                                                       1762, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8938, 3, 1357,
                                                                       1807, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9073, 3, 1393,
                                                                       1852, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9208, 3, 1429,
                                                                       1897, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9343, 3, 1465,
                                                                       1942, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9478, 3, 1537,
                                                                       1987, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9643, 3, 1582,
                                                                       2042, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9808, 3, 1627,
                                                                       2097, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9973, 3, 1672,
                                                                       2152, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10138, 3, 1717,
                                                                       2207, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10303, 3, 1762,
                                                                       2262, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10468, 3, 1807,
                                                                       2317, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10633, 3, 1852,
                                                                       2372, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10798, 3, 1897,
                                                                       2427, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10963, 3, 1987,
                                                                       2482, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11161, 3, 2042,
                                                                       2548, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11359, 3, 2097,
                                                                       2614, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11557, 3, 2152,
                                                                       2680, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11755, 3, 2207,
                                                                       2746, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11953, 3, 2262,
                                                                       2812, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12151, 3, 2317,
                                                                       2878, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12349, 3, 2372,
                                                                       2944, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 12547, 3, 2482,
                                                                       3010, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 12781, 3, 2548,
                                                                       3088, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13015, 3, 2614,
                                                                       3166, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13249, 3, 2680,
                                                                       3244, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13483, 3, 2746,
                                                                       3322, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13717, 3, 2812,
                                                                       3400, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13951, 3, 2878,
                                                                       3478, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14185, 3, 7, 8,
                                                                       3562, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14191, 3, 8, 9,
                                                                       3565, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14197, 3, 9, 10,
                                                                       3568, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14203, 3, 10, 11,
                                                                       3571, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14209, 3, 11, 12,
                                                                       3574, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14215, 3, 12, 13,
                                                                       3577, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14221, 3, 13, 14,
                                                                       3580, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14227, 3, 14, 15,
                                                                       3583, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14233, 3, 15, 16,
                                                                       3586, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14239, 3, 16, 17,
                                                                       3589, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14245, 3, 17, 18,
                                                                       3592, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14251, 3, 18, 19,
                                                                       3595, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14257, 3, 19, 20,
                                                                       3598, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14263, 3, 20, 21,
                                                                       3601, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14269, 3, 21, 22,
                                                                       3604, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14275, 3, 22, 23,
                                                                       3607, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14281, 0, 3,
                                                                       14185, 3562, 14191, 3610,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14299, 0, 3,
                                                                       14191, 3565, 14197, 3619,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14317, 0, 3,
                                                                       14197, 3568, 14203, 3628,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14335, 0, 3,
                                                                       14203, 3571, 14209, 3637,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14353, 0, 3,
                                                                       14209, 3574, 14215, 3646,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14371, 0, 3,
                                                                       14215, 3577, 14221, 3655,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14389, 0, 3,
                                                                       14221, 3580, 14227, 3664,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14407, 0, 3,
                                                                       14227, 3583, 14233, 3673,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14425, 0, 3,
                                                                       14233, 3586, 14239, 3682,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14443, 0, 3,
                                                                       14239, 3589, 14245, 3691,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14461, 0, 3,
                                                                       14245, 3592, 14251, 3700,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14479, 0, 3,
                                                                       14251, 3595, 14257, 3709,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14497, 0, 3,
                                                                       14257, 3598, 14263, 3718,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14515, 0, 3,
                                                                       14263, 3601, 14269, 3727,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 14533, 0, 3,
                                                                       14269, 3604, 14275, 3736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14551, 0, 3,
                                                                       14281, 3610, 14299, 76,
                                                                       82, 3781, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14587, 0, 3,
                                                                       14299, 3619, 14317, 82,
                                                                       88, 3799, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14623, 0, 3,
                                                                       14317, 3628, 14335, 88,
                                                                       94, 3817, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14659, 0, 3,
                                                                       14335, 3637, 14353, 94,
                                                                       100, 3835, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14695, 0, 3,
                                                                       14353, 3646, 14371, 100,
                                                                       106, 3853, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14731, 0, 3,
                                                                       14371, 3655, 14389, 106,
                                                                       112, 3871, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14767, 0, 3,
                                                                       14389, 3664, 14407, 112,
                                                                       118, 3889, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14803, 0, 3,
                                                                       14407, 3673, 14425, 118,
                                                                       124, 3907, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14839, 0, 3,
                                                                       14425, 3682, 14443, 124,
                                                                       130, 3925, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14875, 0, 3,
                                                                       14443, 3691, 14461, 130,
                                                                       136, 3943, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14911, 0, 3,
                                                                       14461, 3700, 14479, 136,
                                                                       142, 3961, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14947, 0, 3,
                                                                       14479, 3709, 14497, 142,
                                                                       148, 3979, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 14983, 0, 3,
                                                                       14497, 3718, 14515, 148,
                                                                       154, 3997, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15019, 0, 3,
                                                                       14515, 3727, 14533, 154,
                                                                       160, 4015, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15055, 0, 3,
                                                                       14551, 3781, 14587, 172,
                                                                       182, 4093, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15115, 0, 3,
                                                                       14587, 3799, 14623, 182,
                                                                       192, 4123, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15175, 0, 3,
                                                                       14623, 3817, 14659, 192,
                                                                       202, 4153, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15235, 0, 3,
                                                                       14659, 3835, 14695, 202,
                                                                       212, 4183, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15295, 0, 3,
                                                                       14695, 3853, 14731, 212,
                                                                       222, 4213, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15355, 0, 3,
                                                                       14731, 3871, 14767, 222,
                                                                       232, 4243, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15415, 0, 3,
                                                                       14767, 3889, 14803, 232,
                                                                       242, 4273, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15475, 0, 3,
                                                                       14803, 3907, 14839, 242,
                                                                       252, 4303, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15535, 0, 3,
                                                                       14839, 3925, 14875, 252,
                                                                       262, 4333, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15595, 0, 3,
                                                                       14875, 3943, 14911, 262,
                                                                       272, 4363, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15655, 0, 3,
                                                                       14911, 3961, 14947, 272,
                                                                       282, 4393, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15715, 0, 3,
                                                                       14947, 3979, 14983, 282,
                                                                       292, 4423, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15775, 0, 3,
                                                                       14983, 3997, 15019, 292,
                                                                       302, 4453, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15835, 0, 3,
                                                                       15055, 4093, 15115, 322,
                                                                       337, 4573, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 15925, 0, 3,
                                                                       15115, 4123, 15175, 337,
                                                                       352, 4618, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16015, 0, 3,
                                                                       15175, 4153, 15235, 352,
                                                                       367, 4663, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16105, 0, 3,
                                                                       15235, 4183, 15295, 367,
                                                                       382, 4708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16195, 0, 3,
                                                                       15295, 4213, 15355, 382,
                                                                       397, 4753, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16285, 0, 3,
                                                                       15355, 4243, 15415, 397,
                                                                       412, 4798, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16375, 0, 3,
                                                                       15415, 4273, 15475, 412,
                                                                       427, 4843, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16465, 0, 3,
                                                                       15475, 4303, 15535, 427,
                                                                       442, 4888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16555, 0, 3,
                                                                       15535, 4333, 15595, 442,
                                                                       457, 4933, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16645, 0, 3,
                                                                       15595, 4363, 15655, 457,
                                                                       472, 4978, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16735, 0, 3,
                                                                       15655, 4393, 15715, 472,
                                                                       487, 5023, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16825, 0, 3,
                                                                       15715, 4423, 15775, 487,
                                                                       502, 5068, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16915, 0, 3,
                                                                       15835, 4573, 15925, 532,
                                                                       553, 5239, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17041, 0, 3,
                                                                       15925, 4618, 16015, 553,
                                                                       574, 5302, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17167, 0, 3,
                                                                       16015, 4663, 16105, 574,
                                                                       595, 5365, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17293, 0, 3,
                                                                       16105, 4708, 16195, 595,
                                                                       616, 5428, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17419, 0, 3,
                                                                       16195, 4753, 16285, 616,
                                                                       637, 5491, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17545, 0, 3,
                                                                       16285, 4798, 16375, 637,
                                                                       658, 5554, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17671, 0, 3,
                                                                       16375, 4843, 16465, 658,
                                                                       679, 5617, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17797, 0, 3,
                                                                       16465, 4888, 16555, 679,
                                                                       700, 5680, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 17923, 0, 3,
                                                                       16555, 4933, 16645, 700,
                                                                       721, 5743, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18049, 0, 3,
                                                                       16645, 4978, 16735, 721,
                                                                       742, 5806, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18175, 0, 3,
                                                                       16735, 5023, 16825, 742,
                                                                       763, 5869, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18301, 0, 3,
                                                                       16915, 5239, 17041, 805,
                                                                       833, 6100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18469, 0, 3,
                                                                       17041, 5302, 17167, 833,
                                                                       861, 6184, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18637, 0, 3,
                                                                       17167, 5365, 17293, 861,
                                                                       889, 6268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18805, 0, 3,
                                                                       17293, 5428, 17419, 889,
                                                                       917, 6352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 18973, 0, 3,
                                                                       17419, 5491, 17545, 917,
                                                                       945, 6436, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19141, 0, 3,
                                                                       17545, 5554, 17671, 945,
                                                                       973, 6520, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19309, 0, 3,
                                                                       17671, 5617, 17797, 973,
                                                                       1001, 6604, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19477, 0, 3,
                                                                       17797, 5680, 17923, 1001,
                                                                       1029, 6688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19645, 0, 3,
                                                                       17923, 5743, 18049, 1029,
                                                                       1057, 6772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19813, 0, 3,
                                                                       18049, 5806, 18175, 1057,
                                                                       1085, 6856, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19981, 0, 3,
                                                                       18301, 6100, 18469, 1141,
                                                                       1177, 7156, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20197, 0, 3,
                                                                       18469, 6184, 18637, 1177,
                                                                       1213, 7264, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20413, 0, 3,
                                                                       18637, 6268, 18805, 1213,
                                                                       1249, 7372, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20629, 0, 3,
                                                                       18805, 6352, 18973, 1249,
                                                                       1285, 7480, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 20845, 0, 3,
                                                                       18973, 6436, 19141, 1285,
                                                                       1321, 7588, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21061, 0, 3,
                                                                       19141, 6520, 19309, 1321,
                                                                       1357, 7696, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21277, 0, 3,
                                                                       19309, 6604, 19477, 1357,
                                                                       1393, 7804, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21493, 0, 3,
                                                                       19477, 6688, 19645, 1393,
                                                                       1429, 7912, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21709, 0, 3,
                                                                       19645, 6772, 19813, 1429,
                                                                       1465, 8020, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21925, 0, 3,
                                                                       19981, 7156, 20197, 1537,
                                                                       1582, 8398, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 22195, 0, 3,
                                                                       20197, 7264, 20413, 1582,
                                                                       1627, 8533, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 22465, 0, 3,
                                                                       20413, 7372, 20629, 1627,
                                                                       1672, 8668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 22735, 0, 3,
                                                                       20629, 7480, 20845, 1672,
                                                                       1717, 8803, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23005, 0, 3,
                                                                       20845, 7588, 21061, 1717,
                                                                       1762, 8938, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23275, 0, 3,
                                                                       21061, 7696, 21277, 1762,
                                                                       1807, 9073, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23545, 0, 3,
                                                                       21277, 7804, 21493, 1807,
                                                                       1852, 9208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 23815, 0, 3,
                                                                       21493, 7912, 21709, 1852,
                                                                       1897, 9343, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 24085, 0, 3,
                                                                       21925, 8398, 22195, 1987,
                                                                       2042, 9808, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 24415, 0, 3,
                                                                       22195, 8533, 22465, 2042,
                                                                       2097, 9973, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 24745, 0, 3,
                                                                       22465, 8668, 22735, 2097,
                                                                       2152, 10138, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 25075, 0, 3,
                                                                       22735, 8803, 23005, 2152,
                                                                       2207, 10303, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 25405, 0, 3,
                                                                       23005, 8938, 23275, 2207,
                                                                       2262, 10468, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 25735, 0, 3,
                                                                       23275, 9073, 23545, 2262,
                                                                       2317, 10633, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 26065, 0, 3,
                                                                       23545, 9208, 23815, 2317,
                                                                       2372, 10798, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 26395, 0, 3,
                                                                       24085, 9808, 24415, 2482,
                                                                       2548, 11359, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 26791, 0, 3,
                                                                       24415, 9973, 24745, 2548,
                                                                       2614, 11557, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 27187, 0, 3,
                                                                       24745, 10138, 25075, 2614,
                                                                       2680, 11755, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 27583, 0, 3,
                                                                       25075, 10303, 25405, 2680,
                                                                       2746, 11953, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 27979, 0, 3,
                                                                       25405, 10468, 25735, 2746,
                                                                       2812, 12151, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 28375, 0, 3,
                                                                       25735, 10633, 26065, 2812,
                                                                       2878, 12349, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 28771, 0, 3,
                                                                       26395, 11359, 26791, 3010,
                                                                       3088, 13015, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 29239, 0, 3,
                                                                       26791, 11557, 27187, 3088,
                                                                       3166, 13249, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 29707, 0, 3,
                                                                       27187, 11755, 27583, 3166,
                                                                       3244, 13483, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 30175, 0, 3,
                                                                       27583, 11953, 27979, 3244,
                                                                       3322, 13717, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 30643, 0, 3,
                                                                       27979, 12151, 28375, 3322,
                                                                       3400, 13951, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31111, 3, 3556,
                                                                       3559, 14185, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31121, 3, 3559,
                                                                       3562, 14191, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31131, 3, 3562,
                                                                       3565, 14197, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31141, 3, 3565,
                                                                       3568, 14203, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31151, 3, 3568,
                                                                       3571, 14209, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31161, 3, 3571,
                                                                       3574, 14215, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31171, 3, 3574,
                                                                       3577, 14221, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31181, 3, 3577,
                                                                       3580, 14227, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31191, 3, 3580,
                                                                       3583, 14233, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31201, 3, 3583,
                                                                       3586, 14239, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31211, 3, 3586,
                                                                       3589, 14245, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31221, 3, 3589,
                                                                       3592, 14251, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31231, 3, 3592,
                                                                       3595, 14257, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31241, 3, 3595,
                                                                       3598, 14263, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31251, 3, 3598,
                                                                       3601, 14269, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 31261, 3, 3601,
                                                                       3604, 14275, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31271, 0, 3,
                                                                       31111, 14185, 31121,
                                                                       14281, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31301, 0, 3,
                                                                       31121, 14191, 31131,
                                                                       14299, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31331, 0, 3,
                                                                       31131, 14197, 31141,
                                                                       14317, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31361, 0, 3,
                                                                       31141, 14203, 31151,
                                                                       14335, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31391, 0, 3,
                                                                       31151, 14209, 31161,
                                                                       14353, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31421, 0, 3,
                                                                       31161, 14215, 31171,
                                                                       14371, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31451, 0, 3,
                                                                       31171, 14221, 31181,
                                                                       14389, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31481, 0, 3,
                                                                       31181, 14227, 31191,
                                                                       14407, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31511, 0, 3,
                                                                       31191, 14233, 31201,
                                                                       14425, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31541, 0, 3,
                                                                       31201, 14239, 31211,
                                                                       14443, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31571, 0, 3,
                                                                       31211, 14245, 31221,
                                                                       14461, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31601, 0, 3,
                                                                       31221, 14251, 31231,
                                                                       14479, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31631, 0, 3,
                                                                       31231, 14257, 31241,
                                                                       14497, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31661, 0, 3,
                                                                       31241, 14263, 31251,
                                                                       14515, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 31691, 0, 3,
                                                                       31251, 14269, 31261,
                                                                       14533, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 31721, 0, 3,
                                                                       31271, 14281, 31301, 3745,
                                                                       3763, 14551, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 31781, 0, 3,
                                                                       31301, 14299, 31331, 3763,
                                                                       3781, 14587, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 31841, 0, 3,
                                                                       31331, 14317, 31361, 3781,
                                                                       3799, 14623, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 31901, 0, 3,
                                                                       31361, 14335, 31391, 3799,
                                                                       3817, 14659, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 31961, 0, 3,
                                                                       31391, 14353, 31421, 3817,
                                                                       3835, 14695, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32021, 0, 3,
                                                                       31421, 14371, 31451, 3835,
                                                                       3853, 14731, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32081, 0, 3,
                                                                       31451, 14389, 31481, 3853,
                                                                       3871, 14767, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32141, 0, 3,
                                                                       31481, 14407, 31511, 3871,
                                                                       3889, 14803, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32201, 0, 3,
                                                                       31511, 14425, 31541, 3889,
                                                                       3907, 14839, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32261, 0, 3,
                                                                       31541, 14443, 31571, 3907,
                                                                       3925, 14875, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32321, 0, 3,
                                                                       31571, 14461, 31601, 3925,
                                                                       3943, 14911, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32381, 0, 3,
                                                                       31601, 14479, 31631, 3943,
                                                                       3961, 14947, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32441, 0, 3,
                                                                       31631, 14497, 31661, 3961,
                                                                       3979, 14983, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 32501, 0, 3,
                                                                       31661, 14515, 31691, 3979,
                                                                       3997, 15019, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 32561, 0, 3,
                                                                       31721, 14551, 31781, 4033,
                                                                       4063, 15055, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 32661, 0, 3,
                                                                       31781, 14587, 31841, 4063,
                                                                       4093, 15115, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 32761, 0, 3,
                                                                       31841, 14623, 31901, 4093,
                                                                       4123, 15175, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 32861, 0, 3,
                                                                       31901, 14659, 31961, 4123,
                                                                       4153, 15235, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 32961, 0, 3,
                                                                       31961, 14695, 32021, 4153,
                                                                       4183, 15295, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33061, 0, 3,
                                                                       32021, 14731, 32081, 4183,
                                                                       4213, 15355, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33161, 0, 3,
                                                                       32081, 14767, 32141, 4213,
                                                                       4243, 15415, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33261, 0, 3,
                                                                       32141, 14803, 32201, 4243,
                                                                       4273, 15475, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33361, 0, 3,
                                                                       32201, 14839, 32261, 4273,
                                                                       4303, 15535, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33461, 0, 3,
                                                                       32261, 14875, 32321, 4303,
                                                                       4333, 15595, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33561, 0, 3,
                                                                       32321, 14911, 32381, 4333,
                                                                       4363, 15655, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33661, 0, 3,
                                                                       32381, 14947, 32441, 4363,
                                                                       4393, 15715, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 33761, 0, 3,
                                                                       32441, 14983, 32501, 4393,
                                                                       4423, 15775, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 33861, 0, 3,
                                                                       32561, 15055, 32661, 4483,
                                                                       4528, 15835, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34011, 0, 3,
                                                                       32661, 15115, 32761, 4528,
                                                                       4573, 15925, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34161, 0, 3,
                                                                       32761, 15175, 32861, 4573,
                                                                       4618, 16015, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34311, 0, 3,
                                                                       32861, 15235, 32961, 4618,
                                                                       4663, 16105, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34461, 0, 3,
                                                                       32961, 15295, 33061, 4663,
                                                                       4708, 16195, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34611, 0, 3,
                                                                       33061, 15355, 33161, 4708,
                                                                       4753, 16285, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34761, 0, 3,
                                                                       33161, 15415, 33261, 4753,
                                                                       4798, 16375, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 34911, 0, 3,
                                                                       33261, 15475, 33361, 4798,
                                                                       4843, 16465, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 35061, 0, 3,
                                                                       33361, 15535, 33461, 4843,
                                                                       4888, 16555, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 35211, 0, 3,
                                                                       33461, 15595, 33561, 4888,
                                                                       4933, 16645, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 35361, 0, 3,
                                                                       33561, 15655, 33661, 4933,
                                                                       4978, 16735, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 35511, 0, 3,
                                                                       33661, 15715, 33761, 4978,
                                                                       5023, 16825, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35661, 0, 3,
                                                                       33861, 15835, 34011, 5113,
                                                                       5176, 16915, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 35871, 0, 3,
                                                                       34011, 15925, 34161, 5176,
                                                                       5239, 17041, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 36081, 0, 3,
                                                                       34161, 16015, 34311, 5239,
                                                                       5302, 17167, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 36291, 0, 3,
                                                                       34311, 16105, 34461, 5302,
                                                                       5365, 17293, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 36501, 0, 3,
                                                                       34461, 16195, 34611, 5365,
                                                                       5428, 17419, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 36711, 0, 3,
                                                                       34611, 16285, 34761, 5428,
                                                                       5491, 17545, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 36921, 0, 3,
                                                                       34761, 16375, 34911, 5491,
                                                                       5554, 17671, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 37131, 0, 3,
                                                                       34911, 16465, 35061, 5554,
                                                                       5617, 17797, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 37341, 0, 3,
                                                                       35061, 16555, 35211, 5617,
                                                                       5680, 17923, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 37551, 0, 3,
                                                                       35211, 16645, 35361, 5680,
                                                                       5743, 18049, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 37761, 0, 3,
                                                                       35361, 16735, 35511, 5743,
                                                                       5806, 18175, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 37971, 0, 3,
                                                                       35661, 16915, 35871, 5932,
                                                                       6016, 18301, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38251, 0, 3,
                                                                       35871, 17041, 36081, 6016,
                                                                       6100, 18469, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38531, 0, 3,
                                                                       36081, 17167, 36291, 6100,
                                                                       6184, 18637, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 38811, 0, 3,
                                                                       36291, 17293, 36501, 6184,
                                                                       6268, 18805, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39091, 0, 3,
                                                                       36501, 17419, 36711, 6268,
                                                                       6352, 18973, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39371, 0, 3,
                                                                       36711, 17545, 36921, 6352,
                                                                       6436, 19141, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39651, 0, 3,
                                                                       36921, 17671, 37131, 6436,
                                                                       6520, 19309, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 39931, 0, 3,
                                                                       37131, 17797, 37341, 6520,
                                                                       6604, 19477, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 40211, 0, 3,
                                                                       37341, 17923, 37551, 6604,
                                                                       6688, 19645, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 40491, 0, 3,
                                                                       37551, 18049, 37761, 6688,
                                                                       6772, 19813, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 40771, 0, 3,
                                                                       37971, 18301, 38251, 6940,
                                                                       7048, 19981, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41131, 0, 3,
                                                                       38251, 18469, 38531, 7048,
                                                                       7156, 20197, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41491, 0, 3,
                                                                       38531, 18637, 38811, 7156,
                                                                       7264, 20413, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 41851, 0, 3,
                                                                       38811, 18805, 39091, 7264,
                                                                       7372, 20629, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42211, 0, 3,
                                                                       39091, 18973, 39371, 7372,
                                                                       7480, 20845, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42571, 0, 3,
                                                                       39371, 19141, 39651, 7480,
                                                                       7588, 21061, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 42931, 0, 3,
                                                                       39651, 19309, 39931, 7588,
                                                                       7696, 21277, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 43291, 0, 3,
                                                                       39931, 19477, 40211, 7696,
                                                                       7804, 21493, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 43651, 0, 3,
                                                                       40211, 19645, 40491, 7804,
                                                                       7912, 21709, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 44011, 0, 3,
                                                                       40771, 19981, 41131, 8128,
                                                                       8263, 21925, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 44461, 0, 3,
                                                                       41131, 20197, 41491, 8263,
                                                                       8398, 22195, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 44911, 0, 3,
                                                                       41491, 20413, 41851, 8398,
                                                                       8533, 22465, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 45361, 0, 3,
                                                                       41851, 20629, 42211, 8533,
                                                                       8668, 22735, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 45811, 0, 3,
                                                                       42211, 20845, 42571, 8668,
                                                                       8803, 23005, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 46261, 0, 3,
                                                                       42571, 21061, 42931, 8803,
                                                                       8938, 23275, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 46711, 0, 3,
                                                                       42931, 21277, 43291, 8938,
                                                                       9073, 23545, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 47161, 0, 3,
                                                                       43291, 21493, 43651, 9073,
                                                                       9208, 23815, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 47611, 0, 3,
                                                                       44011, 21925, 44461, 9478,
                                                                       9643, 24085, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 48161, 0, 3,
                                                                       44461, 22195, 44911, 9643,
                                                                       9808, 24415, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 48711, 0, 3,
                                                                       44911, 22465, 45361, 9808,
                                                                       9973, 24745, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 49261, 0, 3,
                                                                       45361, 22735, 45811, 9973,
                                                                       10138, 25075, ncols,
                                                                       gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 49811, 0, 3,
                                                                       45811, 23005, 46261,
                                                                       10138, 10303, 25405,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 50361, 0, 3,
                                                                       46261, 23275, 46711,
                                                                       10303, 10468, 25735,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 50911, 0, 3,
                                                                       46711, 23545, 47161,
                                                                       10468, 10633, 26065,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 51461, 0, 3,
                                                                       47611, 24085, 48161,
                                                                       10963, 11161, 26395,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 52121, 0, 3,
                                                                       48161, 24415, 48711,
                                                                       11161, 11359, 26791,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 52781, 0, 3,
                                                                       48711, 24745, 49261,
                                                                       11359, 11557, 27187,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 53441, 0, 3,
                                                                       49261, 25075, 49811,
                                                                       11557, 11755, 27583,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 54101, 0, 3,
                                                                       49811, 25405, 50361,
                                                                       11755, 11953, 27979,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 54761, 0, 3,
                                                                       50361, 25735, 50911,
                                                                       11953, 12151, 28375,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 55421, 0, 3,
                                                                       51461, 26395, 52121,
                                                                       12547, 12781, 28771,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 56201, 0, 3,
                                                                       52121, 26791, 52781,
                                                                       12781, 13015, 29239,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 56981, 0, 3,
                                                                       52781, 27187, 53441,
                                                                       13015, 13249, 29707,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 57761, 0, 3,
                                                                       53441, 27583, 54101,
                                                                       13249, 13483, 30175,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 58541, 0, 3,
                                                                       54101, 27979, 54761,
                                                                       13483, 13717, 30643,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59321, 3, 14185,
                                                                       14191, 31131, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59336, 3, 14191,
                                                                       14197, 31141, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59351, 3, 14197,
                                                                       14203, 31151, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59366, 3, 14203,
                                                                       14209, 31161, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59381, 3, 14209,
                                                                       14215, 31171, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59396, 3, 14215,
                                                                       14221, 31181, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59411, 3, 14221,
                                                                       14227, 31191, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59426, 3, 14227,
                                                                       14233, 31201, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59441, 3, 14233,
                                                                       14239, 31211, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59456, 3, 14239,
                                                                       14245, 31221, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59471, 3, 14245,
                                                                       14251, 31231, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59486, 3, 14251,
                                                                       14257, 31241, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59501, 3, 14257,
                                                                       14263, 31251, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 59516, 3, 14263,
                                                                       14269, 31261, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59531, 0, 3,
                                                                       59321, 31131, 59336,
                                                                       14281, 14299, 31331,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59576, 0, 3,
                                                                       59336, 31141, 59351,
                                                                       14299, 14317, 31361,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59621, 0, 3,
                                                                       59351, 31151, 59366,
                                                                       14317, 14335, 31391,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59666, 0, 3,
                                                                       59366, 31161, 59381,
                                                                       14335, 14353, 31421,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59711, 0, 3,
                                                                       59381, 31171, 59396,
                                                                       14353, 14371, 31451,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59756, 0, 3,
                                                                       59396, 31181, 59411,
                                                                       14371, 14389, 31481,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59801, 0, 3,
                                                                       59411, 31191, 59426,
                                                                       14389, 14407, 31511,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59846, 0, 3,
                                                                       59426, 31201, 59441,
                                                                       14407, 14425, 31541,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59891, 0, 3,
                                                                       59441, 31211, 59456,
                                                                       14425, 14443, 31571,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59936, 0, 3,
                                                                       59456, 31221, 59471,
                                                                       14443, 14461, 31601,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 59981, 0, 3,
                                                                       59471, 31231, 59486,
                                                                       14461, 14479, 31631,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 60026, 0, 3,
                                                                       59486, 31241, 59501,
                                                                       14479, 14497, 31661,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 60071, 0, 3,
                                                                       59501, 31251, 59516,
                                                                       14497, 14515, 31691,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60116, 0, 3,
                                                                       59531, 31331, 59576,
                                                                       14551, 14587, 31841,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60206, 0, 3,
                                                                       59576, 31361, 59621,
                                                                       14587, 14623, 31901,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60296, 0, 3,
                                                                       59621, 31391, 59666,
                                                                       14623, 14659, 31961,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60386, 0, 3,
                                                                       59666, 31421, 59711,
                                                                       14659, 14695, 32021,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60476, 0, 3,
                                                                       59711, 31451, 59756,
                                                                       14695, 14731, 32081,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60566, 0, 3,
                                                                       59756, 31481, 59801,
                                                                       14731, 14767, 32141,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60656, 0, 3,
                                                                       59801, 31511, 59846,
                                                                       14767, 14803, 32201,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60746, 0, 3,
                                                                       59846, 31541, 59891,
                                                                       14803, 14839, 32261,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60836, 0, 3,
                                                                       59891, 31571, 59936,
                                                                       14839, 14875, 32321,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 60926, 0, 3,
                                                                       59936, 31601, 59981,
                                                                       14875, 14911, 32381,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 61016, 0, 3,
                                                                       59981, 31631, 60026,
                                                                       14911, 14947, 32441,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 61106, 0, 3,
                                                                       60026, 31661, 60071,
                                                                       14947, 14983, 32501,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61196, 0, 3,
                                                                       60116, 31841, 60206,
                                                                       15055, 15115, 32761,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61346, 0, 3,
                                                                       60206, 31901, 60296,
                                                                       15115, 15175, 32861,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61496, 0, 3,
                                                                       60296, 31961, 60386,
                                                                       15175, 15235, 32961,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61646, 0, 3,
                                                                       60386, 32021, 60476,
                                                                       15235, 15295, 33061,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61796, 0, 3,
                                                                       60476, 32081, 60566,
                                                                       15295, 15355, 33161,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 61946, 0, 3,
                                                                       60566, 32141, 60656,
                                                                       15355, 15415, 33261,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 62096, 0, 3,
                                                                       60656, 32201, 60746,
                                                                       15415, 15475, 33361,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 62246, 0, 3,
                                                                       60746, 32261, 60836,
                                                                       15475, 15535, 33461,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 62396, 0, 3,
                                                                       60836, 32321, 60926,
                                                                       15535, 15595, 33561,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 62546, 0, 3,
                                                                       60926, 32381, 61016,
                                                                       15595, 15655, 33661,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 62696, 0, 3,
                                                                       61016, 32441, 61106,
                                                                       15655, 15715, 33761,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 62846, 0, 3,
                                                                       61196, 32761, 61346,
                                                                       15835, 15925, 34161,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63071, 0, 3,
                                                                       61346, 32861, 61496,
                                                                       15925, 16015, 34311,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63296, 0, 3,
                                                                       61496, 32961, 61646,
                                                                       16015, 16105, 34461,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63521, 0, 3,
                                                                       61646, 33061, 61796,
                                                                       16105, 16195, 34611,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63746, 0, 3,
                                                                       61796, 33161, 61946,
                                                                       16195, 16285, 34761,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 63971, 0, 3,
                                                                       61946, 33261, 62096,
                                                                       16285, 16375, 34911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64196, 0, 3,
                                                                       62096, 33361, 62246,
                                                                       16375, 16465, 35061,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64421, 0, 3,
                                                                       62246, 33461, 62396,
                                                                       16465, 16555, 35211,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64646, 0, 3,
                                                                       62396, 33561, 62546,
                                                                       16555, 16645, 35361,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 64871, 0, 3,
                                                                       62546, 33661, 62696,
                                                                       16645, 16735, 35511,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 65096, 0, 3,
                                                                       62846, 34161, 63071,
                                                                       16915, 17041, 36081,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 65411, 0, 3,
                                                                       63071, 34311, 63296,
                                                                       17041, 17167, 36291,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 65726, 0, 3,
                                                                       63296, 34461, 63521,
                                                                       17167, 17293, 36501,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66041, 0, 3,
                                                                       63521, 34611, 63746,
                                                                       17293, 17419, 36711,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66356, 0, 3,
                                                                       63746, 34761, 63971,
                                                                       17419, 17545, 36921,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66671, 0, 3,
                                                                       63971, 34911, 64196,
                                                                       17545, 17671, 37131,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 66986, 0, 3,
                                                                       64196, 35061, 64421,
                                                                       17671, 17797, 37341,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 67301, 0, 3,
                                                                       64421, 35211, 64646,
                                                                       17797, 17923, 37551,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 67616, 0, 3,
                                                                       64646, 35361, 64871,
                                                                       17923, 18049, 37761,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 67931, 0, 3,
                                                                       65096, 36081, 65411,
                                                                       18301, 18469, 38531,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 68351, 0, 3,
                                                                       65411, 36291, 65726,
                                                                       18469, 18637, 38811,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 68771, 0, 3,
                                                                       65726, 36501, 66041,
                                                                       18637, 18805, 39091,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 69191, 0, 3,
                                                                       66041, 36711, 66356,
                                                                       18805, 18973, 39371,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 69611, 0, 3,
                                                                       66356, 36921, 66671,
                                                                       18973, 19141, 39651,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 70031, 0, 3,
                                                                       66671, 37131, 66986,
                                                                       19141, 19309, 39931,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 70451, 0, 3,
                                                                       66986, 37341, 67301,
                                                                       19309, 19477, 40211,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 70871, 0, 3,
                                                                       67301, 37551, 67616,
                                                                       19477, 19645, 40491,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 71291, 0, 3,
                                                                       67931, 38531, 68351,
                                                                       19981, 20197, 41491,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 71831, 0, 3,
                                                                       68351, 38811, 68771,
                                                                       20197, 20413, 41851,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 72371, 0, 3,
                                                                       68771, 39091, 69191,
                                                                       20413, 20629, 42211,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 72911, 0, 3,
                                                                       69191, 39371, 69611,
                                                                       20629, 20845, 42571,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 73451, 0, 3,
                                                                       69611, 39651, 70031,
                                                                       20845, 21061, 42931,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 73991, 0, 3,
                                                                       70031, 39931, 70451,
                                                                       21061, 21277, 43291,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 74531, 0, 3,
                                                                       70451, 40211, 70871,
                                                                       21277, 21493, 43651,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 75071, 0, 3,
                                                                       71291, 41491, 71831,
                                                                       21925, 22195, 44911,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 75746, 0, 3,
                                                                       71831, 41851, 72371,
                                                                       22195, 22465, 45361,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 76421, 0, 3,
                                                                       72371, 42211, 72911,
                                                                       22465, 22735, 45811,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 77096, 0, 3,
                                                                       72911, 42571, 73451,
                                                                       22735, 23005, 46261,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 77771, 0, 3,
                                                                       73451, 42931, 73991,
                                                                       23005, 23275, 46711,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 78446, 0, 3,
                                                                       73991, 43291, 74531,
                                                                       23275, 23545, 47161,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 79121, 0, 3,
                                                                       75071, 44911, 75746,
                                                                       24085, 24415, 48711,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 79946, 0, 3,
                                                                       75746, 45361, 76421,
                                                                       24415, 24745, 49261,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 80771, 0, 3,
                                                                       76421, 45811, 77096,
                                                                       24745, 25075, 49811,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 81596, 0, 3,
                                                                       77096, 46261, 77771,
                                                                       25075, 25405, 50361,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 82421, 0, 3,
                                                                       77771, 46711, 78446,
                                                                       25405, 25735, 50911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 83246, 0, 3,
                                                                       79121, 48711, 79946,
                                                                       26395, 26791, 52781,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 84236, 0, 3,
                                                                       79946, 49261, 80771,
                                                                       26791, 27187, 53441,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 85226, 0, 3,
                                                                       80771, 49811, 81596,
                                                                       27187, 27583, 54101,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 86216, 0, 3,
                                                                       81596, 50361, 82421,
                                                                       27583, 27979, 54761,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 87206, 0, 3,
                                                                       83246, 52781, 84236,
                                                                       28771, 29239, 56981,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 88376, 0, 3,
                                                                       84236, 53441, 85226,
                                                                       29239, 29707, 57761,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 89546, 0, 3,
                                                                       85226, 54101, 86216,
                                                                       29707, 30175, 58541,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90716, 3, 31111,
                                                                       31121, 59321, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90737, 3, 31121,
                                                                       31131, 59336, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90758, 3, 31131,
                                                                       31141, 59351, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90779, 3, 31141,
                                                                       31151, 59366, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90800, 3, 31151,
                                                                       31161, 59381, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90821, 3, 31161,
                                                                       31171, 59396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90842, 3, 31171,
                                                                       31181, 59411, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90863, 3, 31181,
                                                                       31191, 59426, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90884, 3, 31191,
                                                                       31201, 59441, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90905, 3, 31201,
                                                                       31211, 59456, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90926, 3, 31211,
                                                                       31221, 59471, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90947, 3, 31221,
                                                                       31231, 59486, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90968, 3, 31231,
                                                                       31241, 59501, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 90989, 3, 31241,
                                                                       31251, 59516, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91010, 0, 3,
                                                                       90716, 59321, 90737,
                                                                       31271, 31301, 59531,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91073, 0, 3,
                                                                       90737, 59336, 90758,
                                                                       31301, 31331, 59576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91136, 0, 3,
                                                                       90758, 59351, 90779,
                                                                       31331, 31361, 59621,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91199, 0, 3,
                                                                       90779, 59366, 90800,
                                                                       31361, 31391, 59666,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91262, 0, 3,
                                                                       90800, 59381, 90821,
                                                                       31391, 31421, 59711,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91325, 0, 3,
                                                                       90821, 59396, 90842,
                                                                       31421, 31451, 59756,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91388, 0, 3,
                                                                       90842, 59411, 90863,
                                                                       31451, 31481, 59801,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91451, 0, 3,
                                                                       90863, 59426, 90884,
                                                                       31481, 31511, 59846,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91514, 0, 3,
                                                                       90884, 59441, 90905,
                                                                       31511, 31541, 59891,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91577, 0, 3,
                                                                       90905, 59456, 90926,
                                                                       31541, 31571, 59936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91640, 0, 3,
                                                                       90926, 59471, 90947,
                                                                       31571, 31601, 59981,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91703, 0, 3,
                                                                       90947, 59486, 90968,
                                                                       31601, 31631, 60026,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 91766, 0, 3,
                                                                       90968, 59501, 90989,
                                                                       31631, 31661, 60071,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 91829, 0, 3,
                                                                       91010, 59531, 91073,
                                                                       31721, 31781, 60116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 91955, 0, 3,
                                                                       91073, 59576, 91136,
                                                                       31781, 31841, 60206,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92081, 0, 3,
                                                                       91136, 59621, 91199,
                                                                       31841, 31901, 60296,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92207, 0, 3,
                                                                       91199, 59666, 91262,
                                                                       31901, 31961, 60386,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92333, 0, 3,
                                                                       91262, 59711, 91325,
                                                                       31961, 32021, 60476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92459, 0, 3,
                                                                       91325, 59756, 91388,
                                                                       32021, 32081, 60566,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92585, 0, 3,
                                                                       91388, 59801, 91451,
                                                                       32081, 32141, 60656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92711, 0, 3,
                                                                       91451, 59846, 91514,
                                                                       32141, 32201, 60746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92837, 0, 3,
                                                                       91514, 59891, 91577,
                                                                       32201, 32261, 60836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 92963, 0, 3,
                                                                       91577, 59936, 91640,
                                                                       32261, 32321, 60926,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93089, 0, 3,
                                                                       91640, 59981, 91703,
                                                                       32321, 32381, 61016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 93215, 0, 3,
                                                                       91703, 60026, 91766,
                                                                       32381, 32441, 61106,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 93341, 0, 3,
                                                                       91829, 60116, 91955,
                                                                       32561, 32661, 61196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 93551, 0, 3,
                                                                       91955, 60206, 92081,
                                                                       32661, 32761, 61346,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 93761, 0, 3,
                                                                       92081, 60296, 92207,
                                                                       32761, 32861, 61496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 93971, 0, 3,
                                                                       92207, 60386, 92333,
                                                                       32861, 32961, 61646,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94181, 0, 3,
                                                                       92333, 60476, 92459,
                                                                       32961, 33061, 61796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94391, 0, 3,
                                                                       92459, 60566, 92585,
                                                                       33061, 33161, 61946,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94601, 0, 3,
                                                                       92585, 60656, 92711,
                                                                       33161, 33261, 62096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 94811, 0, 3,
                                                                       92711, 60746, 92837,
                                                                       33261, 33361, 62246,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95021, 0, 3,
                                                                       92837, 60836, 92963,
                                                                       33361, 33461, 62396,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95231, 0, 3,
                                                                       92963, 60926, 93089,
                                                                       33461, 33561, 62546,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 95441, 0, 3,
                                                                       93089, 61016, 93215,
                                                                       33561, 33661, 62696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 95651, 0, 3,
                                                                       93341, 61196, 93551,
                                                                       33861, 34011, 62846,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 95966, 0, 3,
                                                                       93551, 61346, 93761,
                                                                       34011, 34161, 63071,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 96281, 0, 3,
                                                                       93761, 61496, 93971,
                                                                       34161, 34311, 63296,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 96596, 0, 3,
                                                                       93971, 61646, 94181,
                                                                       34311, 34461, 63521,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 96911, 0, 3,
                                                                       94181, 61796, 94391,
                                                                       34461, 34611, 63746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 97226, 0, 3,
                                                                       94391, 61946, 94601,
                                                                       34611, 34761, 63971,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 97541, 0, 3,
                                                                       94601, 62096, 94811,
                                                                       34761, 34911, 64196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 97856, 0, 3,
                                                                       94811, 62246, 95021,
                                                                       34911, 35061, 64421,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 98171, 0, 3,
                                                                       95021, 62396, 95231,
                                                                       35061, 35211, 64646,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 98486, 0, 3,
                                                                       95231, 62546, 95441,
                                                                       35211, 35361, 64871,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 98801, 0, 3,
                                                                       95651, 62846, 95966,
                                                                       35661, 35871, 65096,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 99242, 0, 3,
                                                                       95966, 63071, 96281,
                                                                       35871, 36081, 65411,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 99683, 0, 3,
                                                                       96281, 63296, 96596,
                                                                       36081, 36291, 65726,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 100124, 0, 3,
                                                                       96596, 63521, 96911,
                                                                       36291, 36501, 66041,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 100565, 0, 3,
                                                                       96911, 63746, 97226,
                                                                       36501, 36711, 66356,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 101006, 0, 3,
                                                                       97226, 63971, 97541,
                                                                       36711, 36921, 66671,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 101447, 0, 3,
                                                                       97541, 64196, 97856,
                                                                       36921, 37131, 66986,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 101888, 0, 3,
                                                                       97856, 64421, 98171,
                                                                       37131, 37341, 67301,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 102329, 0, 3,
                                                                       98171, 64646, 98486,
                                                                       37341, 37551, 67616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 102770, 0, 3,
                                                                       98801, 65096, 99242,
                                                                       37971, 38251, 67931,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 103358, 0, 3,
                                                                       99242, 65411, 99683,
                                                                       38251, 38531, 68351,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 103946, 0, 3,
                                                                       99683, 65726, 100124,
                                                                       38531, 38811, 68771,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 104534, 0, 3,
                                                                       100124, 66041, 100565,
                                                                       38811, 39091, 69191,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 105122, 0, 3,
                                                                       100565, 66356, 101006,
                                                                       39091, 39371, 69611,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 105710, 0, 3,
                                                                       101006, 66671, 101447,
                                                                       39371, 39651, 70031,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 106298, 0, 3,
                                                                       101447, 66986, 101888,
                                                                       39651, 39931, 70451,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 106886, 0, 3,
                                                                       101888, 67301, 102329,
                                                                       39931, 40211, 70871,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 107474, 0, 3,
                                                                       102770, 67931, 103358,
                                                                       40771, 41131, 71291,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 108230, 0, 3,
                                                                       103358, 68351, 103946,
                                                                       41131, 41491, 71831,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 108986, 0, 3,
                                                                       103946, 68771, 104534,
                                                                       41491, 41851, 72371,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 109742, 0, 3,
                                                                       104534, 69191, 105122,
                                                                       41851, 42211, 72911,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 110498, 0, 3,
                                                                       105122, 69611, 105710,
                                                                       42211, 42571, 73451,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 111254, 0, 3,
                                                                       105710, 70031, 106298,
                                                                       42571, 42931, 73991,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 112010, 0, 3,
                                                                       106298, 70451, 106886,
                                                                       42931, 43291, 74531,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 112766, 0, 3,
                                                                       107474, 71291, 108230,
                                                                       44011, 44461, 75071,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 113711, 0, 3,
                                                                       108230, 71831, 108986,
                                                                       44461, 44911, 75746,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 114656, 0, 3,
                                                                       108986, 72371, 109742,
                                                                       44911, 45361, 76421,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 115601, 0, 3,
                                                                       109742, 72911, 110498,
                                                                       45361, 45811, 77096,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 116546, 0, 3,
                                                                       110498, 73451, 111254,
                                                                       45811, 46261, 77771,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 117491, 0, 3,
                                                                       111254, 73991, 112010,
                                                                       46261, 46711, 78446,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 118436, 0, 3,
                                                                       112766, 75071, 113711,
                                                                       47611, 48161, 79121,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 119591, 0, 3,
                                                                       113711, 75746, 114656,
                                                                       48161, 48711, 79946,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 120746, 0, 3,
                                                                       114656, 76421, 115601,
                                                                       48711, 49261, 80771,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 121901, 0, 3,
                                                                       115601, 77096, 116546,
                                                                       49261, 49811, 81596,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 123056, 0, 3,
                                                                       116546, 77771, 117491,
                                                                       49811, 50361, 82421,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 124211, 0, 3,
                                                                       118436, 79121, 119591,
                                                                       51461, 52121, 83246,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 125597, 0, 3,
                                                                       119591, 79946, 120746,
                                                                       52121, 52781, 84236,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 126983, 0, 3,
                                                                       120746, 80771, 121901,
                                                                       52781, 53441, 85226,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 128369, 0, 3,
                                                                       121901, 81596, 123056,
                                                                       53441, 54101, 86216,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 129755, 0, 3,
                                                                       124211, 83246, 125597,
                                                                       55421, 56201, 87206,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 131393, 0, 3,
                                                                       125597, 84236, 126983,
                                                                       56201, 56981, 88376,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 133031, 0, 3,
                                                                       126983, 85226, 128369,
                                                                       56981, 57761, 89546,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134669, 3, 59321,
                                                                       59336, 90758, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134697, 3, 59336,
                                                                       59351, 90779, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134725, 3, 59351,
                                                                       59366, 90800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134753, 3, 59366,
                                                                       59381, 90821, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134781, 3, 59381,
                                                                       59396, 90842, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134809, 3, 59396,
                                                                       59411, 90863, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134837, 3, 59411,
                                                                       59426, 90884, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134865, 3, 59426,
                                                                       59441, 90905, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134893, 3, 59441,
                                                                       59456, 90926, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134921, 3, 59456,
                                                                       59471, 90947, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134949, 3, 59471,
                                                                       59486, 90968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 134977, 3, 59486,
                                                                       59501, 90989, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135005, 0, 3,
                                                                       134669, 90758, 134697,
                                                                       59531, 59576, 91136,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135089, 0, 3,
                                                                       134697, 90779, 134725,
                                                                       59576, 59621, 91199,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135173, 0, 3,
                                                                       134725, 90800, 134753,
                                                                       59621, 59666, 91262,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135257, 0, 3,
                                                                       134753, 90821, 134781,
                                                                       59666, 59711, 91325,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135341, 0, 3,
                                                                       134781, 90842, 134809,
                                                                       59711, 59756, 91388,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135425, 0, 3,
                                                                       134809, 90863, 134837,
                                                                       59756, 59801, 91451,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135509, 0, 3,
                                                                       134837, 90884, 134865,
                                                                       59801, 59846, 91514,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135593, 0, 3,
                                                                       134865, 90905, 134893,
                                                                       59846, 59891, 91577,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135677, 0, 3,
                                                                       134893, 90926, 134921,
                                                                       59891, 59936, 91640,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135761, 0, 3,
                                                                       134921, 90947, 134949,
                                                                       59936, 59981, 91703,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 135845, 0, 3,
                                                                       134949, 90968, 134977,
                                                                       59981, 60026, 91766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 135929, 0, 3,
                                                                       135005, 91136, 135089,
                                                                       60116, 60206, 92081,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136097, 0, 3,
                                                                       135089, 91199, 135173,
                                                                       60206, 60296, 92207,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136265, 0, 3,
                                                                       135173, 91262, 135257,
                                                                       60296, 60386, 92333,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136433, 0, 3,
                                                                       135257, 91325, 135341,
                                                                       60386, 60476, 92459,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136601, 0, 3,
                                                                       135341, 91388, 135425,
                                                                       60476, 60566, 92585,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136769, 0, 3,
                                                                       135425, 91451, 135509,
                                                                       60566, 60656, 92711,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 136937, 0, 3,
                                                                       135509, 91514, 135593,
                                                                       60656, 60746, 92837,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 137105, 0, 3,
                                                                       135593, 91577, 135677,
                                                                       60746, 60836, 92963,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 137273, 0, 3,
                                                                       135677, 91640, 135761,
                                                                       60836, 60926, 93089,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 137441, 0, 3,
                                                                       135761, 91703, 135845,
                                                                       60926, 61016, 93215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 137609, 0, 3,
                                                                       135929, 92081, 136097,
                                                                       61196, 61346, 93761,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 137889, 0, 3,
                                                                       136097, 92207, 136265,
                                                                       61346, 61496, 93971,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 138169, 0, 3,
                                                                       136265, 92333, 136433,
                                                                       61496, 61646, 94181,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 138449, 0, 3,
                                                                       136433, 92459, 136601,
                                                                       61646, 61796, 94391,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 138729, 0, 3,
                                                                       136601, 92585, 136769,
                                                                       61796, 61946, 94601,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 139009, 0, 3,
                                                                       136769, 92711, 136937,
                                                                       61946, 62096, 94811,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 139289, 0, 3,
                                                                       136937, 92837, 137105,
                                                                       62096, 62246, 95021,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 139569, 0, 3,
                                                                       137105, 92963, 137273,
                                                                       62246, 62396, 95231,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 139849, 0, 3,
                                                                       137273, 93089, 137441,
                                                                       62396, 62546, 95441,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 140129, 0, 3,
                                                                       137609, 93761, 137889,
                                                                       62846, 63071, 96281,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 140549, 0, 3,
                                                                       137889, 93971, 138169,
                                                                       63071, 63296, 96596,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 140969, 0, 3,
                                                                       138169, 94181, 138449,
                                                                       63296, 63521, 96911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 141389, 0, 3,
                                                                       138449, 94391, 138729,
                                                                       63521, 63746, 97226,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 141809, 0, 3,
                                                                       138729, 94601, 139009,
                                                                       63746, 63971, 97541,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 142229, 0, 3,
                                                                       139009, 94811, 139289,
                                                                       63971, 64196, 97856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 142649, 0, 3,
                                                                       139289, 95021, 139569,
                                                                       64196, 64421, 98171,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 143069, 0, 3,
                                                                       139569, 95231, 139849,
                                                                       64421, 64646, 98486,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 143489, 0, 3,
                                                                       140129, 96281, 140549,
                                                                       65096, 65411, 99683,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 144077, 0, 3,
                                                                       140549, 96596, 140969,
                                                                       65411, 65726, 100124,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 144665, 0, 3,
                                                                       140969, 96911, 141389,
                                                                       65726, 66041, 100565,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 145253, 0, 3,
                                                                       141389, 97226, 141809,
                                                                       66041, 66356, 101006,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 145841, 0, 3,
                                                                       141809, 97541, 142229,
                                                                       66356, 66671, 101447,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 146429, 0, 3,
                                                                       142229, 97856, 142649,
                                                                       66671, 66986, 101888,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 147017, 0, 3,
                                                                       142649, 98171, 143069,
                                                                       66986, 67301, 102329,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 147605, 0, 3,
                                                                       143489, 99683, 144077,
                                                                       67931, 68351, 103946,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 148389, 0, 3,
                                                                       144077, 100124, 144665,
                                                                       68351, 68771, 104534,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 149173, 0, 3,
                                                                       144665, 100565, 145253,
                                                                       68771, 69191, 105122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 149957, 0, 3,
                                                                       145253, 101006, 145841,
                                                                       69191, 69611, 105710,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 150741, 0, 3,
                                                                       145841, 101447, 146429,
                                                                       69611, 70031, 106298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 151525, 0, 3,
                                                                       146429, 101888, 147017,
                                                                       70031, 70451, 106886,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 152309, 0, 3,
                                                                       147605, 103946, 148389,
                                                                       71291, 71831, 108986,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 153317, 0, 3,
                                                                       148389, 104534, 149173,
                                                                       71831, 72371, 109742,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 154325, 0, 3,
                                                                       149173, 105122, 149957,
                                                                       72371, 72911, 110498,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 155333, 0, 3,
                                                                       149957, 105710, 150741,
                                                                       72911, 73451, 111254,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 156341, 0, 3,
                                                                       150741, 106298, 151525,
                                                                       73451, 73991, 112010,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 157349, 0, 3,
                                                                       152309, 108986, 153317,
                                                                       75071, 75746, 114656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 158609, 0, 3,
                                                                       153317, 109742, 154325,
                                                                       75746, 76421, 115601,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 159869, 0, 3,
                                                                       154325, 110498, 155333,
                                                                       76421, 77096, 116546,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 161129, 0, 3,
                                                                       155333, 111254, 156341,
                                                                       77096, 77771, 117491,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 162389, 0, 3,
                                                                       157349, 114656, 158609,
                                                                       79121, 79946, 120746,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 163929, 0, 3,
                                                                       158609, 115601, 159869,
                                                                       79946, 80771, 121901,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 165469, 0, 3,
                                                                       159869, 116546, 161129,
                                                                       80771, 81596, 123056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 167009, 0, 3,
                                                                       162389, 120746, 163929,
                                                                       83246, 84236, 126983,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 168857, 0, 3,
                                                                       163929, 121901, 165469,
                                                                       84236, 85226, 128369,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 170705, 0, 3,
                                                                       167009, 126983, 168857,
                                                                       87206, 88376, 133031,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172889, 3, 90716,
                                                                       90737, 134669, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172925, 3, 90737,
                                                                       90758, 134697, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172961, 3, 90758,
                                                                       90779, 134725, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 172997, 3, 90779,
                                                                       90800, 134753, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173033, 3, 90800,
                                                                       90821, 134781, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173069, 3, 90821,
                                                                       90842, 134809, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173105, 3, 90842,
                                                                       90863, 134837, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173141, 3, 90863,
                                                                       90884, 134865, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173177, 3, 90884,
                                                                       90905, 134893, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173213, 3, 90905,
                                                                       90926, 134921, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173249, 3, 90926,
                                                                       90947, 134949, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 173285, 3, 90947,
                                                                       90968, 134977, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173321, 0, 3,
                                                                       172889, 134669, 172925,
                                                                       91010, 91073, 135005,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173429, 0, 3,
                                                                       172925, 134697, 172961,
                                                                       91073, 91136, 135089,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173537, 0, 3,
                                                                       172961, 134725, 172997,
                                                                       91136, 91199, 135173,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173645, 0, 3,
                                                                       172997, 134753, 173033,
                                                                       91199, 91262, 135257,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173753, 0, 3,
                                                                       173033, 134781, 173069,
                                                                       91262, 91325, 135341,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173861, 0, 3,
                                                                       173069, 134809, 173105,
                                                                       91325, 91388, 135425,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 173969, 0, 3,
                                                                       173105, 134837, 173141,
                                                                       91388, 91451, 135509,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 174077, 0, 3,
                                                                       173141, 134865, 173177,
                                                                       91451, 91514, 135593,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 174185, 0, 3,
                                                                       173177, 134893, 173213,
                                                                       91514, 91577, 135677,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 174293, 0, 3,
                                                                       173213, 134921, 173249,
                                                                       91577, 91640, 135761,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 174401, 0, 3,
                                                                       173249, 134949, 173285,
                                                                       91640, 91703, 135845,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 174509, 0, 3,
                                                                       173321, 135005, 173429,
                                                                       91829, 91955, 135929,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 174725, 0, 3,
                                                                       173429, 135089, 173537,
                                                                       91955, 92081, 136097,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 174941, 0, 3,
                                                                       173537, 135173, 173645,
                                                                       92081, 92207, 136265,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 175157, 0, 3,
                                                                       173645, 135257, 173753,
                                                                       92207, 92333, 136433,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 175373, 0, 3,
                                                                       173753, 135341, 173861,
                                                                       92333, 92459, 136601,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 175589, 0, 3,
                                                                       173861, 135425, 173969,
                                                                       92459, 92585, 136769,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 175805, 0, 3,
                                                                       173969, 135509, 174077,
                                                                       92585, 92711, 136937,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 176021, 0, 3,
                                                                       174077, 135593, 174185,
                                                                       92711, 92837, 137105,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 176237, 0, 3,
                                                                       174185, 135677, 174293,
                                                                       92837, 92963, 137273,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 176453, 0, 3,
                                                                       174293, 135761, 174401,
                                                                       92963, 93089, 137441,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 176669, 0, 3,
                                                                       174509, 135929, 174725,
                                                                       93341, 93551, 137609,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 177029, 0, 3,
                                                                       174725, 136097, 174941,
                                                                       93551, 93761, 137889,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 177389, 0, 3,
                                                                       174941, 136265, 175157,
                                                                       93761, 93971, 138169,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 177749, 0, 3,
                                                                       175157, 136433, 175373,
                                                                       93971, 94181, 138449,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 178109, 0, 3,
                                                                       175373, 136601, 175589,
                                                                       94181, 94391, 138729,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 178469, 0, 3,
                                                                       175589, 136769, 175805,
                                                                       94391, 94601, 139009,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 178829, 0, 3,
                                                                       175805, 136937, 176021,
                                                                       94601, 94811, 139289,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 179189, 0, 3,
                                                                       176021, 137105, 176237,
                                                                       94811, 95021, 139569,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 179549, 0, 3,
                                                                       176237, 137273, 176453,
                                                                       95021, 95231, 139849,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 179909, 0, 3,
                                                                       176669, 137609, 177029,
                                                                       95651, 95966, 140129,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 180449, 0, 3,
                                                                       177029, 137889, 177389,
                                                                       95966, 96281, 140549,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 180989, 0, 3,
                                                                       177389, 138169, 177749,
                                                                       96281, 96596, 140969,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 181529, 0, 3,
                                                                       177749, 138449, 178109,
                                                                       96596, 96911, 141389,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 182069, 0, 3,
                                                                       178109, 138729, 178469,
                                                                       96911, 97226, 141809,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 182609, 0, 3,
                                                                       178469, 139009, 178829,
                                                                       97226, 97541, 142229,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 183149, 0, 3,
                                                                       178829, 139289, 179189,
                                                                       97541, 97856, 142649,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 183689, 0, 3,
                                                                       179189, 139569, 179549,
                                                                       97856, 98171, 143069,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 184229, 0, 3,
                                                                       179909, 140129, 180449,
                                                                       98801, 99242, 143489,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 184985, 0, 3,
                                                                       180449, 140549, 180989,
                                                                       99242, 99683, 144077,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 185741, 0, 3,
                                                                       180989, 140969, 181529,
                                                                       99683, 100124, 144665,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 186497, 0, 3,
                                                                       181529, 141389, 182069,
                                                                       100124, 100565, 145253,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 187253, 0, 3,
                                                                       182069, 141809, 182609,
                                                                       100565, 101006, 145841,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 188009, 0, 3,
                                                                       182609, 142229, 183149,
                                                                       101006, 101447, 146429,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 188765, 0, 3,
                                                                       183149, 142649, 183689,
                                                                       101447, 101888, 147017,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 189521, 0, 3,
                                                                       184229, 143489, 184985,
                                                                       102770, 103358, 147605,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 190529, 0, 3,
                                                                       184985, 144077, 185741,
                                                                       103358, 103946, 148389,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 191537, 0, 3,
                                                                       185741, 144665, 186497,
                                                                       103946, 104534, 149173,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 192545, 0, 3,
                                                                       186497, 145253, 187253,
                                                                       104534, 105122, 149957,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 193553, 0, 3,
                                                                       187253, 145841, 188009,
                                                                       105122, 105710, 150741,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 194561, 0, 3,
                                                                       188009, 146429, 188765,
                                                                       105710, 106298, 151525,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 195569, 0, 3,
                                                                       189521, 147605, 190529,
                                                                       107474, 108230, 152309,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 196865, 0, 3,
                                                                       190529, 148389, 191537,
                                                                       108230, 108986, 153317,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 198161, 0, 3,
                                                                       191537, 149173, 192545,
                                                                       108986, 109742, 154325,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 199457, 0, 3,
                                                                       192545, 149957, 193553,
                                                                       109742, 110498, 155333,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 200753, 0, 3,
                                                                       193553, 150741, 194561,
                                                                       110498, 111254, 156341,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 202049, 0, 3,
                                                                       195569, 152309, 196865,
                                                                       112766, 113711, 157349,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 203669, 0, 3,
                                                                       196865, 153317, 198161,
                                                                       113711, 114656, 158609,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 205289, 0, 3,
                                                                       198161, 154325, 199457,
                                                                       114656, 115601, 159869,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 206909, 0, 3,
                                                                       199457, 155333, 200753,
                                                                       115601, 116546, 161129,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 208529, 0, 3,
                                                                       202049, 157349, 203669,
                                                                       118436, 119591, 162389,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 210509, 0, 3,
                                                                       203669, 158609, 205289,
                                                                       119591, 120746, 163929,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 212489, 0, 3,
                                                                       205289, 159869, 206909,
                                                                       120746, 121901, 165469,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 214469, 0, 3,
                                                                       208529, 162389, 210509,
                                                                       124211, 125597, 167009,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 216845, 0, 3,
                                                                       210509, 163929, 212489,
                                                                       125597, 126983, 168857,
                                                                       ncols, gamma, p, q);

                    compute_prim_sok_three_center_electron_repulsion_0(buffer, 219221, 0, 3,
                                                                       214469, 167009, 216845,
                                                                       129755, 131393, 170705,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 222029, 189521, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 223457, 195569, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 225293, 202049, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 227588, 208529, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 230393, 214469, 2376, ncols);

                    simdfunc::contract_primitives(buffer, 233759, 219221, 2808, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 223037, 222029, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 224753, 223457, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 226913, 225293, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 229568, 227588, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 232769, 230393, 66, 1, nmax);

        simdtrf::transform_k_inner(buffer, 236567, 233759, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 237737, 223037, 224753, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 238997, 224753, 226913, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 240617, 226913, 229568, 15, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 242642, 229568, 232769, 15, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 245117, 232769, 236567, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 248087, 237737, 238997, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 250607, 238997, 240617, 15, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 253847, 240617, 242642, 15, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 257897, 242642, 245117, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 262847, 248087, 250607, 15, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 267047, 250607, 253847, 15, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 272447, 253847, 257897, 15, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 279197, 262847, 267047, 15, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 285497, 267047, 272447, 15, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 293597, 279197, 285497, 15, nmax);

        simdtrf::transform_i_inner(buffer, 302417, 293597, 21, 15, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 302417, 195, nmax);
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
