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


#include "SimdThreeCenterElectronRepulsionRecHII.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
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
compute_hii_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hii_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 220388, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1859 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 220388, 148155, 11614, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 17,
                                                             ncols, fj, mu, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3556, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3559, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3562, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3565, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3568, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3571, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3574, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3577, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3580, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3583, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3586, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3589, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3592, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3595, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3598, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 3601, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3604, 3, 9, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3613, 3, 10, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3622, 3, 11, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3631, 3, 12, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3640, 3, 13, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3649, 3, 14, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3658, 3, 15, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3667, 3, 16, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3676, 3, 17, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3685, 3, 18, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3694, 3, 19, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3703, 3, 20, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3712, 3, 21, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3721, 3, 22, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 3730, 3, 23, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3739, 3, 31, 88,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3757, 3, 34, 94,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3775, 3, 37, 100,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3793, 3, 40, 106,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3811, 3, 43, 112,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3829, 3, 46, 118,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3847, 3, 49, 124,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3865, 3, 52, 130,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3883, 3, 55, 136,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3901, 3, 58, 142,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3919, 3, 61, 148,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3937, 3, 64, 154,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3955, 3, 67, 160,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3973, 3, 70, 166,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3991, 3, 88, 192,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4021, 3, 94, 202,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4051, 3, 100, 212,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4081, 3, 106, 222,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4111, 3, 112, 232,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4141, 3, 118, 242,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4171, 3, 124, 252,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4201, 3, 130, 262,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4231, 3, 136, 272,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4261, 3, 142, 282,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4291, 3, 148, 292,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4321, 3, 154, 302,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4351, 3, 160, 312,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4381, 3, 192, 352,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4426, 3, 202, 367,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4471, 3, 212, 382,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4516, 3, 222, 397,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4561, 3, 232, 412,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4606, 3, 242, 427,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4651, 3, 252, 442,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4696, 3, 262, 457,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4741, 3, 272, 472,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4786, 3, 282, 487,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4831, 3, 292, 502,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4876, 3, 302, 517,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4921, 3, 352, 574,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4984, 3, 367, 595,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5047, 3, 382, 616,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5110, 3, 397, 637,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5173, 3, 412, 658,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5236, 3, 427, 679,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5299, 3, 442, 700,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5362, 3, 457, 721,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5425, 3, 472, 742,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5488, 3, 487, 763,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5551, 3, 502, 784,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5614, 3, 574, 861,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5698, 3, 595, 889,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5782, 3, 616, 917,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5866, 3, 637, 945,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5950, 3, 658, 973,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6034, 3, 679,
                                                                       1001, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6118, 3, 700,
                                                                       1029, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6202, 3, 721,
                                                                       1057, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6286, 3, 742,
                                                                       1085, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6370, 3, 763,
                                                                       1113, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6454, 3, 861,
                                                                       1213, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6562, 3, 889,
                                                                       1249, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6670, 3, 917,
                                                                       1285, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6778, 3, 945,
                                                                       1321, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6886, 3, 973,
                                                                       1357, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6994, 3, 1001,
                                                                       1393, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7102, 3, 1029,
                                                                       1429, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7210, 3, 1057,
                                                                       1465, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7318, 3, 1085,
                                                                       1501, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7426, 3, 1213,
                                                                       1627, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7561, 3, 1249,
                                                                       1672, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7696, 3, 1285,
                                                                       1717, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7831, 3, 1321,
                                                                       1762, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7966, 3, 1357,
                                                                       1807, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8101, 3, 1393,
                                                                       1852, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8236, 3, 1429,
                                                                       1897, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8371, 3, 1465,
                                                                       1942, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8506, 3, 1627,
                                                                       2097, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8671, 3, 1672,
                                                                       2152, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8836, 3, 1717,
                                                                       2207, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9001, 3, 1762,
                                                                       2262, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9166, 3, 1807,
                                                                       2317, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9331, 3, 1852,
                                                                       2372, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9496, 3, 1897,
                                                                       2427, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9661, 3, 2097,
                                                                       2614, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9859, 3, 2152,
                                                                       2680, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10057, 3, 2207,
                                                                       2746, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10255, 3, 2262,
                                                                       2812, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10453, 3, 2317,
                                                                       2878, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10651, 3, 2372,
                                                                       2944, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 10849, 3, 2614,
                                                                       3166, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 11083, 3, 2680,
                                                                       3244, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 11317, 3, 2746,
                                                                       3322, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 11551, 3, 2812,
                                                                       3400, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 11785, 3, 2878,
                                                                       3478, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12019, 3, 7, 8,
                                                                       3556, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12025, 3, 8, 9,
                                                                       3559, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12031, 3, 9, 10,
                                                                       3562, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12037, 3, 10, 11,
                                                                       3565, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12043, 3, 11, 12,
                                                                       3568, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12049, 3, 12, 13,
                                                                       3571, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12055, 3, 13, 14,
                                                                       3574, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12061, 3, 14, 15,
                                                                       3577, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12067, 3, 15, 16,
                                                                       3580, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12073, 3, 16, 17,
                                                                       3583, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12079, 3, 17, 18,
                                                                       3586, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12085, 3, 18, 19,
                                                                       3589, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12091, 3, 19, 20,
                                                                       3592, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12097, 3, 20, 21,
                                                                       3595, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12103, 3, 21, 22,
                                                                       3598, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 12109, 3, 22, 23,
                                                                       3601, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12115, 0, 3,
                                                                       12019, 3556, 12025, 3604,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12133, 0, 3,
                                                                       12025, 3559, 12031, 3613,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12151, 0, 3,
                                                                       12031, 3562, 12037, 3622,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12169, 0, 3,
                                                                       12037, 3565, 12043, 3631,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12187, 0, 3,
                                                                       12043, 3568, 12049, 3640,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12205, 0, 3,
                                                                       12049, 3571, 12055, 3649,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12223, 0, 3,
                                                                       12055, 3574, 12061, 3658,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12241, 0, 3,
                                                                       12061, 3577, 12067, 3667,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12259, 0, 3,
                                                                       12067, 3580, 12073, 3676,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12277, 0, 3,
                                                                       12073, 3583, 12079, 3685,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12295, 0, 3,
                                                                       12079, 3586, 12085, 3694,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12313, 0, 3,
                                                                       12085, 3589, 12091, 3703,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12331, 0, 3,
                                                                       12091, 3592, 12097, 3712,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12349, 0, 3,
                                                                       12097, 3595, 12103, 3721,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 12367, 0, 3,
                                                                       12103, 3598, 12109, 3730,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12385, 0, 3,
                                                                       12115, 3604, 12133, 76,
                                                                       82, 3739, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12421, 0, 3,
                                                                       12133, 3613, 12151, 82,
                                                                       88, 3757, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12457, 0, 3,
                                                                       12151, 3622, 12169, 88,
                                                                       94, 3775, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12493, 0, 3,
                                                                       12169, 3631, 12187, 94,
                                                                       100, 3793, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12529, 0, 3,
                                                                       12187, 3640, 12205, 100,
                                                                       106, 3811, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12565, 0, 3,
                                                                       12205, 3649, 12223, 106,
                                                                       112, 3829, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12601, 0, 3,
                                                                       12223, 3658, 12241, 112,
                                                                       118, 3847, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12637, 0, 3,
                                                                       12241, 3667, 12259, 118,
                                                                       124, 3865, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12673, 0, 3,
                                                                       12259, 3676, 12277, 124,
                                                                       130, 3883, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12709, 0, 3,
                                                                       12277, 3685, 12295, 130,
                                                                       136, 3901, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12745, 0, 3,
                                                                       12295, 3694, 12313, 136,
                                                                       142, 3919, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12781, 0, 3,
                                                                       12313, 3703, 12331, 142,
                                                                       148, 3937, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12817, 0, 3,
                                                                       12331, 3712, 12349, 148,
                                                                       154, 3955, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 12853, 0, 3,
                                                                       12349, 3721, 12367, 154,
                                                                       160, 3973, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12889, 0, 3,
                                                                       12385, 3739, 12421, 172,
                                                                       182, 3991, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12949, 0, 3,
                                                                       12421, 3757, 12457, 182,
                                                                       192, 4021, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13009, 0, 3,
                                                                       12457, 3775, 12493, 192,
                                                                       202, 4051, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13069, 0, 3,
                                                                       12493, 3793, 12529, 202,
                                                                       212, 4081, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13129, 0, 3,
                                                                       12529, 3811, 12565, 212,
                                                                       222, 4111, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13189, 0, 3,
                                                                       12565, 3829, 12601, 222,
                                                                       232, 4141, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13249, 0, 3,
                                                                       12601, 3847, 12637, 232,
                                                                       242, 4171, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13309, 0, 3,
                                                                       12637, 3865, 12673, 242,
                                                                       252, 4201, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13369, 0, 3,
                                                                       12673, 3883, 12709, 252,
                                                                       262, 4231, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13429, 0, 3,
                                                                       12709, 3901, 12745, 262,
                                                                       272, 4261, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13489, 0, 3,
                                                                       12745, 3919, 12781, 272,
                                                                       282, 4291, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13549, 0, 3,
                                                                       12781, 3937, 12817, 282,
                                                                       292, 4321, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 13609, 0, 3,
                                                                       12817, 3955, 12853, 292,
                                                                       302, 4351, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13669, 0, 3,
                                                                       12889, 3991, 12949, 322,
                                                                       337, 4381, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13759, 0, 3,
                                                                       12949, 4021, 13009, 337,
                                                                       352, 4426, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13849, 0, 3,
                                                                       13009, 4051, 13069, 352,
                                                                       367, 4471, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13939, 0, 3,
                                                                       13069, 4081, 13129, 367,
                                                                       382, 4516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14029, 0, 3,
                                                                       13129, 4111, 13189, 382,
                                                                       397, 4561, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14119, 0, 3,
                                                                       13189, 4141, 13249, 397,
                                                                       412, 4606, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14209, 0, 3,
                                                                       13249, 4171, 13309, 412,
                                                                       427, 4651, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14299, 0, 3,
                                                                       13309, 4201, 13369, 427,
                                                                       442, 4696, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14389, 0, 3,
                                                                       13369, 4231, 13429, 442,
                                                                       457, 4741, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14479, 0, 3,
                                                                       13429, 4261, 13489, 457,
                                                                       472, 4786, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14569, 0, 3,
                                                                       13489, 4291, 13549, 472,
                                                                       487, 4831, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 14659, 0, 3,
                                                                       13549, 4321, 13609, 487,
                                                                       502, 4876, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14749, 0, 3,
                                                                       13669, 4381, 13759, 532,
                                                                       553, 4921, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14875, 0, 3,
                                                                       13759, 4426, 13849, 553,
                                                                       574, 4984, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15001, 0, 3,
                                                                       13849, 4471, 13939, 574,
                                                                       595, 5047, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15127, 0, 3,
                                                                       13939, 4516, 14029, 595,
                                                                       616, 5110, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15253, 0, 3,
                                                                       14029, 4561, 14119, 616,
                                                                       637, 5173, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15379, 0, 3,
                                                                       14119, 4606, 14209, 637,
                                                                       658, 5236, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15505, 0, 3,
                                                                       14209, 4651, 14299, 658,
                                                                       679, 5299, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15631, 0, 3,
                                                                       14299, 4696, 14389, 679,
                                                                       700, 5362, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15757, 0, 3,
                                                                       14389, 4741, 14479, 700,
                                                                       721, 5425, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 15883, 0, 3,
                                                                       14479, 4786, 14569, 721,
                                                                       742, 5488, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 16009, 0, 3,
                                                                       14569, 4831, 14659, 742,
                                                                       763, 5551, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16135, 0, 3,
                                                                       14749, 4921, 14875, 805,
                                                                       833, 5614, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16303, 0, 3,
                                                                       14875, 4984, 15001, 833,
                                                                       861, 5698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16471, 0, 3,
                                                                       15001, 5047, 15127, 861,
                                                                       889, 5782, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16639, 0, 3,
                                                                       15127, 5110, 15253, 889,
                                                                       917, 5866, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16807, 0, 3,
                                                                       15253, 5173, 15379, 917,
                                                                       945, 5950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 16975, 0, 3,
                                                                       15379, 5236, 15505, 945,
                                                                       973, 6034, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17143, 0, 3,
                                                                       15505, 5299, 15631, 973,
                                                                       1001, 6118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17311, 0, 3,
                                                                       15631, 5362, 15757, 1001,
                                                                       1029, 6202, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17479, 0, 3,
                                                                       15757, 5425, 15883, 1029,
                                                                       1057, 6286, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 17647, 0, 3,
                                                                       15883, 5488, 16009, 1057,
                                                                       1085, 6370, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17815, 0, 3,
                                                                       16135, 5614, 16303, 1141,
                                                                       1177, 6454, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18031, 0, 3,
                                                                       16303, 5698, 16471, 1177,
                                                                       1213, 6562, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18247, 0, 3,
                                                                       16471, 5782, 16639, 1213,
                                                                       1249, 6670, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18463, 0, 3,
                                                                       16639, 5866, 16807, 1249,
                                                                       1285, 6778, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18679, 0, 3,
                                                                       16807, 5950, 16975, 1285,
                                                                       1321, 6886, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 18895, 0, 3,
                                                                       16975, 6034, 17143, 1321,
                                                                       1357, 6994, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19111, 0, 3,
                                                                       17143, 6118, 17311, 1357,
                                                                       1393, 7102, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19327, 0, 3,
                                                                       17311, 6202, 17479, 1393,
                                                                       1429, 7210, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 19543, 0, 3,
                                                                       17479, 6286, 17647, 1429,
                                                                       1465, 7318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 19759, 0, 3,
                                                                       17815, 6454, 18031, 1537,
                                                                       1582, 7426, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 20029, 0, 3,
                                                                       18031, 6562, 18247, 1582,
                                                                       1627, 7561, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 20299, 0, 3,
                                                                       18247, 6670, 18463, 1627,
                                                                       1672, 7696, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 20569, 0, 3,
                                                                       18463, 6778, 18679, 1672,
                                                                       1717, 7831, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 20839, 0, 3,
                                                                       18679, 6886, 18895, 1717,
                                                                       1762, 7966, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21109, 0, 3,
                                                                       18895, 6994, 19111, 1762,
                                                                       1807, 8101, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21379, 0, 3,
                                                                       19111, 7102, 19327, 1807,
                                                                       1852, 8236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 21649, 0, 3,
                                                                       19327, 7210, 19543, 1852,
                                                                       1897, 8371, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 21919, 0, 3,
                                                                       19759, 7426, 20029, 1987,
                                                                       2042, 8506, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 22249, 0, 3,
                                                                       20029, 7561, 20299, 2042,
                                                                       2097, 8671, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 22579, 0, 3,
                                                                       20299, 7696, 20569, 2097,
                                                                       2152, 8836, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 22909, 0, 3,
                                                                       20569, 7831, 20839, 2152,
                                                                       2207, 9001, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 23239, 0, 3,
                                                                       20839, 7966, 21109, 2207,
                                                                       2262, 9166, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 23569, 0, 3,
                                                                       21109, 8101, 21379, 2262,
                                                                       2317, 9331, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 23899, 0, 3,
                                                                       21379, 8236, 21649, 2317,
                                                                       2372, 9496, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 24229, 0, 3,
                                                                       21919, 8506, 22249, 2482,
                                                                       2548, 9661, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 24625, 0, 3,
                                                                       22249, 8671, 22579, 2548,
                                                                       2614, 9859, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 25021, 0, 3,
                                                                       22579, 8836, 22909, 2614,
                                                                       2680, 10057, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 25417, 0, 3,
                                                                       22909, 9001, 23239, 2680,
                                                                       2746, 10255, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 25813, 0, 3,
                                                                       23239, 9166, 23569, 2746,
                                                                       2812, 10453, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 26209, 0, 3,
                                                                       23569, 9331, 23899, 2812,
                                                                       2878, 10651, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 26605, 0, 3,
                                                                       24229, 9661, 24625, 3010,
                                                                       3088, 10849, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 27073, 0, 3,
                                                                       24625, 9859, 25021, 3088,
                                                                       3166, 11083, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 27541, 0, 3,
                                                                       25021, 10057, 25417, 3166,
                                                                       3244, 11317, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 28009, 0, 3,
                                                                       25417, 10255, 25813, 3244,
                                                                       3322, 11551, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 28477, 0, 3,
                                                                       25813, 10453, 26209, 3322,
                                                                       3400, 11785, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28945, 3, 3556,
                                                                       3559, 12031, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28955, 3, 3559,
                                                                       3562, 12037, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28965, 3, 3562,
                                                                       3565, 12043, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28975, 3, 3565,
                                                                       3568, 12049, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28985, 3, 3568,
                                                                       3571, 12055, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 28995, 3, 3571,
                                                                       3574, 12061, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29005, 3, 3574,
                                                                       3577, 12067, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29015, 3, 3577,
                                                                       3580, 12073, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29025, 3, 3580,
                                                                       3583, 12079, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29035, 3, 3583,
                                                                       3586, 12085, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29045, 3, 3586,
                                                                       3589, 12091, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29055, 3, 3589,
                                                                       3592, 12097, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29065, 3, 3592,
                                                                       3595, 12103, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 29075, 3, 3595,
                                                                       3598, 12109, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29085, 0, 3,
                                                                       28945, 12031, 28955,
                                                                       12151, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29115, 0, 3,
                                                                       28955, 12037, 28965,
                                                                       12169, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29145, 0, 3,
                                                                       28965, 12043, 28975,
                                                                       12187, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29175, 0, 3,
                                                                       28975, 12049, 28985,
                                                                       12205, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29205, 0, 3,
                                                                       28985, 12055, 28995,
                                                                       12223, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29235, 0, 3,
                                                                       28995, 12061, 29005,
                                                                       12241, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29265, 0, 3,
                                                                       29005, 12067, 29015,
                                                                       12259, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29295, 0, 3,
                                                                       29015, 12073, 29025,
                                                                       12277, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29325, 0, 3,
                                                                       29025, 12079, 29035,
                                                                       12295, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29355, 0, 3,
                                                                       29035, 12085, 29045,
                                                                       12313, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29385, 0, 3,
                                                                       29045, 12091, 29055,
                                                                       12331, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29415, 0, 3,
                                                                       29055, 12097, 29065,
                                                                       12349, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 29445, 0, 3,
                                                                       29065, 12103, 29075,
                                                                       12367, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 29475, 0, 3,
                                                                       29085, 12151, 29115, 3739,
                                                                       3757, 12457, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 29535, 0, 3,
                                                                       29115, 12169, 29145, 3757,
                                                                       3775, 12493, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 29595, 0, 3,
                                                                       29145, 12187, 29175, 3775,
                                                                       3793, 12529, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 29655, 0, 3,
                                                                       29175, 12205, 29205, 3793,
                                                                       3811, 12565, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 29715, 0, 3,
                                                                       29205, 12223, 29235, 3811,
                                                                       3829, 12601, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 29775, 0, 3,
                                                                       29235, 12241, 29265, 3829,
                                                                       3847, 12637, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 29835, 0, 3,
                                                                       29265, 12259, 29295, 3847,
                                                                       3865, 12673, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 29895, 0, 3,
                                                                       29295, 12277, 29325, 3865,
                                                                       3883, 12709, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 29955, 0, 3,
                                                                       29325, 12295, 29355, 3883,
                                                                       3901, 12745, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 30015, 0, 3,
                                                                       29355, 12313, 29385, 3901,
                                                                       3919, 12781, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 30075, 0, 3,
                                                                       29385, 12331, 29415, 3919,
                                                                       3937, 12817, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 30135, 0, 3,
                                                                       29415, 12349, 29445, 3937,
                                                                       3955, 12853, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30195, 0, 3,
                                                                       29475, 12457, 29535, 3991,
                                                                       4021, 13009, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30295, 0, 3,
                                                                       29535, 12493, 29595, 4021,
                                                                       4051, 13069, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30395, 0, 3,
                                                                       29595, 12529, 29655, 4051,
                                                                       4081, 13129, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30495, 0, 3,
                                                                       29655, 12565, 29715, 4081,
                                                                       4111, 13189, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30595, 0, 3,
                                                                       29715, 12601, 29775, 4111,
                                                                       4141, 13249, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30695, 0, 3,
                                                                       29775, 12637, 29835, 4141,
                                                                       4171, 13309, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30795, 0, 3,
                                                                       29835, 12673, 29895, 4171,
                                                                       4201, 13369, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30895, 0, 3,
                                                                       29895, 12709, 29955, 4201,
                                                                       4231, 13429, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 30995, 0, 3,
                                                                       29955, 12745, 30015, 4231,
                                                                       4261, 13489, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 31095, 0, 3,
                                                                       30015, 12781, 30075, 4261,
                                                                       4291, 13549, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 31195, 0, 3,
                                                                       30075, 12817, 30135, 4291,
                                                                       4321, 13609, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31295, 0, 3,
                                                                       30195, 13009, 30295, 4381,
                                                                       4426, 13849, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31445, 0, 3,
                                                                       30295, 13069, 30395, 4426,
                                                                       4471, 13939, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31595, 0, 3,
                                                                       30395, 13129, 30495, 4471,
                                                                       4516, 14029, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31745, 0, 3,
                                                                       30495, 13189, 30595, 4516,
                                                                       4561, 14119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 31895, 0, 3,
                                                                       30595, 13249, 30695, 4561,
                                                                       4606, 14209, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32045, 0, 3,
                                                                       30695, 13309, 30795, 4606,
                                                                       4651, 14299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32195, 0, 3,
                                                                       30795, 13369, 30895, 4651,
                                                                       4696, 14389, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32345, 0, 3,
                                                                       30895, 13429, 30995, 4696,
                                                                       4741, 14479, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32495, 0, 3,
                                                                       30995, 13489, 31095, 4741,
                                                                       4786, 14569, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 32645, 0, 3,
                                                                       31095, 13549, 31195, 4786,
                                                                       4831, 14659, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 32795, 0, 3,
                                                                       31295, 13849, 31445, 4921,
                                                                       4984, 15001, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33005, 0, 3,
                                                                       31445, 13939, 31595, 4984,
                                                                       5047, 15127, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33215, 0, 3,
                                                                       31595, 14029, 31745, 5047,
                                                                       5110, 15253, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33425, 0, 3,
                                                                       31745, 14119, 31895, 5110,
                                                                       5173, 15379, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33635, 0, 3,
                                                                       31895, 14209, 32045, 5173,
                                                                       5236, 15505, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 33845, 0, 3,
                                                                       32045, 14299, 32195, 5236,
                                                                       5299, 15631, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34055, 0, 3,
                                                                       32195, 14389, 32345, 5299,
                                                                       5362, 15757, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34265, 0, 3,
                                                                       32345, 14479, 32495, 5362,
                                                                       5425, 15883, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 34475, 0, 3,
                                                                       32495, 14569, 32645, 5425,
                                                                       5488, 16009, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 34685, 0, 3,
                                                                       32795, 15001, 33005, 5614,
                                                                       5698, 16471, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 34965, 0, 3,
                                                                       33005, 15127, 33215, 5698,
                                                                       5782, 16639, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 35245, 0, 3,
                                                                       33215, 15253, 33425, 5782,
                                                                       5866, 16807, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 35525, 0, 3,
                                                                       33425, 15379, 33635, 5866,
                                                                       5950, 16975, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 35805, 0, 3,
                                                                       33635, 15505, 33845, 5950,
                                                                       6034, 17143, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36085, 0, 3,
                                                                       33845, 15631, 34055, 6034,
                                                                       6118, 17311, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36365, 0, 3,
                                                                       34055, 15757, 34265, 6118,
                                                                       6202, 17479, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 36645, 0, 3,
                                                                       34265, 15883, 34475, 6202,
                                                                       6286, 17647, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 36925, 0, 3,
                                                                       34685, 16471, 34965, 6454,
                                                                       6562, 18247, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 37285, 0, 3,
                                                                       34965, 16639, 35245, 6562,
                                                                       6670, 18463, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 37645, 0, 3,
                                                                       35245, 16807, 35525, 6670,
                                                                       6778, 18679, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 38005, 0, 3,
                                                                       35525, 16975, 35805, 6778,
                                                                       6886, 18895, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 38365, 0, 3,
                                                                       35805, 17143, 36085, 6886,
                                                                       6994, 19111, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 38725, 0, 3,
                                                                       36085, 17311, 36365, 6994,
                                                                       7102, 19327, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 39085, 0, 3,
                                                                       36365, 17479, 36645, 7102,
                                                                       7210, 19543, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 39445, 0, 3,
                                                                       36925, 18247, 37285, 7426,
                                                                       7561, 20299, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 39895, 0, 3,
                                                                       37285, 18463, 37645, 7561,
                                                                       7696, 20569, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 40345, 0, 3,
                                                                       37645, 18679, 38005, 7696,
                                                                       7831, 20839, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 40795, 0, 3,
                                                                       38005, 18895, 38365, 7831,
                                                                       7966, 21109, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 41245, 0, 3,
                                                                       38365, 19111, 38725, 7966,
                                                                       8101, 21379, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 41695, 0, 3,
                                                                       38725, 19327, 39085, 8101,
                                                                       8236, 21649, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 42145, 0, 3,
                                                                       39445, 20299, 39895, 8506,
                                                                       8671, 22579, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 42695, 0, 3,
                                                                       39895, 20569, 40345, 8671,
                                                                       8836, 22909, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 43245, 0, 3,
                                                                       40345, 20839, 40795, 8836,
                                                                       9001, 23239, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 43795, 0, 3,
                                                                       40795, 21109, 41245, 9001,
                                                                       9166, 23569, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 44345, 0, 3,
                                                                       41245, 21379, 41695, 9166,
                                                                       9331, 23899, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 44895, 0, 3,
                                                                       42145, 22579, 42695, 9661,
                                                                       9859, 25021, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 45555, 0, 3,
                                                                       42695, 22909, 43245, 9859,
                                                                       10057, 25417, ncols,
                                                                       gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 46215, 0, 3,
                                                                       43245, 23239, 43795,
                                                                       10057, 10255, 25813,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 46875, 0, 3,
                                                                       43795, 23569, 44345,
                                                                       10255, 10453, 26209,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 47535, 0, 3,
                                                                       44895, 25021, 45555,
                                                                       10849, 11083, 27541,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 48315, 0, 3,
                                                                       45555, 25417, 46215,
                                                                       11083, 11317, 28009,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 49095, 0, 3,
                                                                       46215, 25813, 46875,
                                                                       11317, 11551, 28477,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49875, 3, 12019,
                                                                       12025, 28945, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49890, 3, 12025,
                                                                       12031, 28955, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49905, 3, 12031,
                                                                       12037, 28965, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49920, 3, 12037,
                                                                       12043, 28975, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49935, 3, 12043,
                                                                       12049, 28985, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49950, 3, 12049,
                                                                       12055, 28995, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49965, 3, 12055,
                                                                       12061, 29005, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49980, 3, 12061,
                                                                       12067, 29015, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 49995, 3, 12067,
                                                                       12073, 29025, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50010, 3, 12073,
                                                                       12079, 29035, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50025, 3, 12079,
                                                                       12085, 29045, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50040, 3, 12085,
                                                                       12091, 29055, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50055, 3, 12091,
                                                                       12097, 29065, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 50070, 3, 12097,
                                                                       12103, 29075, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50085, 0, 3,
                                                                       49875, 28945, 49890,
                                                                       12115, 12133, 29085,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50130, 0, 3,
                                                                       49890, 28955, 49905,
                                                                       12133, 12151, 29115,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50175, 0, 3,
                                                                       49905, 28965, 49920,
                                                                       12151, 12169, 29145,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50220, 0, 3,
                                                                       49920, 28975, 49935,
                                                                       12169, 12187, 29175,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50265, 0, 3,
                                                                       49935, 28985, 49950,
                                                                       12187, 12205, 29205,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50310, 0, 3,
                                                                       49950, 28995, 49965,
                                                                       12205, 12223, 29235,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50355, 0, 3,
                                                                       49965, 29005, 49980,
                                                                       12223, 12241, 29265,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50400, 0, 3,
                                                                       49980, 29015, 49995,
                                                                       12241, 12259, 29295,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50445, 0, 3,
                                                                       49995, 29025, 50010,
                                                                       12259, 12277, 29325,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50490, 0, 3,
                                                                       50010, 29035, 50025,
                                                                       12277, 12295, 29355,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50535, 0, 3,
                                                                       50025, 29045, 50040,
                                                                       12295, 12313, 29385,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50580, 0, 3,
                                                                       50040, 29055, 50055,
                                                                       12313, 12331, 29415,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 50625, 0, 3,
                                                                       50055, 29065, 50070,
                                                                       12331, 12349, 29445,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50670, 0, 3,
                                                                       50085, 29085, 50130,
                                                                       12385, 12421, 29475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50760, 0, 3,
                                                                       50130, 29115, 50175,
                                                                       12421, 12457, 29535,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50850, 0, 3,
                                                                       50175, 29145, 50220,
                                                                       12457, 12493, 29595,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 50940, 0, 3,
                                                                       50220, 29175, 50265,
                                                                       12493, 12529, 29655,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51030, 0, 3,
                                                                       50265, 29205, 50310,
                                                                       12529, 12565, 29715,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51120, 0, 3,
                                                                       50310, 29235, 50355,
                                                                       12565, 12601, 29775,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51210, 0, 3,
                                                                       50355, 29265, 50400,
                                                                       12601, 12637, 29835,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51300, 0, 3,
                                                                       50400, 29295, 50445,
                                                                       12637, 12673, 29895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51390, 0, 3,
                                                                       50445, 29325, 50490,
                                                                       12673, 12709, 29955,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51480, 0, 3,
                                                                       50490, 29355, 50535,
                                                                       12709, 12745, 30015,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51570, 0, 3,
                                                                       50535, 29385, 50580,
                                                                       12745, 12781, 30075,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 51660, 0, 3,
                                                                       50580, 29415, 50625,
                                                                       12781, 12817, 30135,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 51750, 0, 3,
                                                                       50670, 29475, 50760,
                                                                       12889, 12949, 30195,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 51900, 0, 3,
                                                                       50760, 29535, 50850,
                                                                       12949, 13009, 30295,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52050, 0, 3,
                                                                       50850, 29595, 50940,
                                                                       13009, 13069, 30395,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52200, 0, 3,
                                                                       50940, 29655, 51030,
                                                                       13069, 13129, 30495,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52350, 0, 3,
                                                                       51030, 29715, 51120,
                                                                       13129, 13189, 30595,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52500, 0, 3,
                                                                       51120, 29775, 51210,
                                                                       13189, 13249, 30695,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52650, 0, 3,
                                                                       51210, 29835, 51300,
                                                                       13249, 13309, 30795,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52800, 0, 3,
                                                                       51300, 29895, 51390,
                                                                       13309, 13369, 30895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 52950, 0, 3,
                                                                       51390, 29955, 51480,
                                                                       13369, 13429, 30995,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53100, 0, 3,
                                                                       51480, 30015, 51570,
                                                                       13429, 13489, 31095,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 53250, 0, 3,
                                                                       51570, 30075, 51660,
                                                                       13489, 13549, 31195,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 53400, 0, 3,
                                                                       51750, 30195, 51900,
                                                                       13669, 13759, 31295,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 53625, 0, 3,
                                                                       51900, 30295, 52050,
                                                                       13759, 13849, 31445,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 53850, 0, 3,
                                                                       52050, 30395, 52200,
                                                                       13849, 13939, 31595,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54075, 0, 3,
                                                                       52200, 30495, 52350,
                                                                       13939, 14029, 31745,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54300, 0, 3,
                                                                       52350, 30595, 52500,
                                                                       14029, 14119, 31895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54525, 0, 3,
                                                                       52500, 30695, 52650,
                                                                       14119, 14209, 32045,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54750, 0, 3,
                                                                       52650, 30795, 52800,
                                                                       14209, 14299, 32195,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 54975, 0, 3,
                                                                       52800, 30895, 52950,
                                                                       14299, 14389, 32345,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55200, 0, 3,
                                                                       52950, 30995, 53100,
                                                                       14389, 14479, 32495,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 55425, 0, 3,
                                                                       53100, 31095, 53250,
                                                                       14479, 14569, 32645,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 55650, 0, 3,
                                                                       53400, 31295, 53625,
                                                                       14749, 14875, 32795,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 55965, 0, 3,
                                                                       53625, 31445, 53850,
                                                                       14875, 15001, 33005,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 56280, 0, 3,
                                                                       53850, 31595, 54075,
                                                                       15001, 15127, 33215,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 56595, 0, 3,
                                                                       54075, 31745, 54300,
                                                                       15127, 15253, 33425,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 56910, 0, 3,
                                                                       54300, 31895, 54525,
                                                                       15253, 15379, 33635,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 57225, 0, 3,
                                                                       54525, 32045, 54750,
                                                                       15379, 15505, 33845,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 57540, 0, 3,
                                                                       54750, 32195, 54975,
                                                                       15505, 15631, 34055,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 57855, 0, 3,
                                                                       54975, 32345, 55200,
                                                                       15631, 15757, 34265,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 58170, 0, 3,
                                                                       55200, 32495, 55425,
                                                                       15757, 15883, 34475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 58485, 0, 3,
                                                                       55650, 32795, 55965,
                                                                       16135, 16303, 34685,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 58905, 0, 3,
                                                                       55965, 33005, 56280,
                                                                       16303, 16471, 34965,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 59325, 0, 3,
                                                                       56280, 33215, 56595,
                                                                       16471, 16639, 35245,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 59745, 0, 3,
                                                                       56595, 33425, 56910,
                                                                       16639, 16807, 35525,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 60165, 0, 3,
                                                                       56910, 33635, 57225,
                                                                       16807, 16975, 35805,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 60585, 0, 3,
                                                                       57225, 33845, 57540,
                                                                       16975, 17143, 36085,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 61005, 0, 3,
                                                                       57540, 34055, 57855,
                                                                       17143, 17311, 36365,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 61425, 0, 3,
                                                                       57855, 34265, 58170,
                                                                       17311, 17479, 36645,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 61845, 0, 3,
                                                                       58485, 34685, 58905,
                                                                       17815, 18031, 36925,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 62385, 0, 3,
                                                                       58905, 34965, 59325,
                                                                       18031, 18247, 37285,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 62925, 0, 3,
                                                                       59325, 35245, 59745,
                                                                       18247, 18463, 37645,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 63465, 0, 3,
                                                                       59745, 35525, 60165,
                                                                       18463, 18679, 38005,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 64005, 0, 3,
                                                                       60165, 35805, 60585,
                                                                       18679, 18895, 38365,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 64545, 0, 3,
                                                                       60585, 36085, 61005,
                                                                       18895, 19111, 38725,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 65085, 0, 3,
                                                                       61005, 36365, 61425,
                                                                       19111, 19327, 39085,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 65625, 0, 3,
                                                                       61845, 36925, 62385,
                                                                       19759, 20029, 39445,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 66300, 0, 3,
                                                                       62385, 37285, 62925,
                                                                       20029, 20299, 39895,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 66975, 0, 3,
                                                                       62925, 37645, 63465,
                                                                       20299, 20569, 40345,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 67650, 0, 3,
                                                                       63465, 38005, 64005,
                                                                       20569, 20839, 40795,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 68325, 0, 3,
                                                                       64005, 38365, 64545,
                                                                       20839, 21109, 41245,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 69000, 0, 3,
                                                                       64545, 38725, 65085,
                                                                       21109, 21379, 41695,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 69675, 0, 3,
                                                                       65625, 39445, 66300,
                                                                       21919, 22249, 42145,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 70500, 0, 3,
                                                                       66300, 39895, 66975,
                                                                       22249, 22579, 42695,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 71325, 0, 3,
                                                                       66975, 40345, 67650,
                                                                       22579, 22909, 43245,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 72150, 0, 3,
                                                                       67650, 40795, 68325,
                                                                       22909, 23239, 43795,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 72975, 0, 3,
                                                                       68325, 41245, 69000,
                                                                       23239, 23569, 44345,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 73800, 0, 3,
                                                                       69675, 42145, 70500,
                                                                       24229, 24625, 44895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 74790, 0, 3,
                                                                       70500, 42695, 71325,
                                                                       24625, 25021, 45555,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 75780, 0, 3,
                                                                       71325, 43245, 72150,
                                                                       25021, 25417, 46215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 76770, 0, 3,
                                                                       72150, 43795, 72975,
                                                                       25417, 25813, 46875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 77760, 0, 3,
                                                                       73800, 44895, 74790,
                                                                       26605, 27073, 47535,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 78930, 0, 3,
                                                                       74790, 45555, 75780,
                                                                       27073, 27541, 48315,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 80100, 0, 3,
                                                                       75780, 46215, 76770,
                                                                       27541, 28009, 49095,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81270, 3, 28945,
                                                                       28955, 49905, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81291, 3, 28955,
                                                                       28965, 49920, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81312, 3, 28965,
                                                                       28975, 49935, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81333, 3, 28975,
                                                                       28985, 49950, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81354, 3, 28985,
                                                                       28995, 49965, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81375, 3, 28995,
                                                                       29005, 49980, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81396, 3, 29005,
                                                                       29015, 49995, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81417, 3, 29015,
                                                                       29025, 50010, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81438, 3, 29025,
                                                                       29035, 50025, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81459, 3, 29035,
                                                                       29045, 50040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81480, 3, 29045,
                                                                       29055, 50055, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 81501, 3, 29055,
                                                                       29065, 50070, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 81522, 0, 3,
                                                                       81270, 49905, 81291,
                                                                       29085, 29115, 50175,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 81585, 0, 3,
                                                                       81291, 49920, 81312,
                                                                       29115, 29145, 50220,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 81648, 0, 3,
                                                                       81312, 49935, 81333,
                                                                       29145, 29175, 50265,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 81711, 0, 3,
                                                                       81333, 49950, 81354,
                                                                       29175, 29205, 50310,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 81774, 0, 3,
                                                                       81354, 49965, 81375,
                                                                       29205, 29235, 50355,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 81837, 0, 3,
                                                                       81375, 49980, 81396,
                                                                       29235, 29265, 50400,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 81900, 0, 3,
                                                                       81396, 49995, 81417,
                                                                       29265, 29295, 50445,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 81963, 0, 3,
                                                                       81417, 50010, 81438,
                                                                       29295, 29325, 50490,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 82026, 0, 3,
                                                                       81438, 50025, 81459,
                                                                       29325, 29355, 50535,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 82089, 0, 3,
                                                                       81459, 50040, 81480,
                                                                       29355, 29385, 50580,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 82152, 0, 3,
                                                                       81480, 50055, 81501,
                                                                       29385, 29415, 50625,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 82215, 0, 3,
                                                                       81522, 50175, 81585,
                                                                       29475, 29535, 50850,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 82341, 0, 3,
                                                                       81585, 50220, 81648,
                                                                       29535, 29595, 50940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 82467, 0, 3,
                                                                       81648, 50265, 81711,
                                                                       29595, 29655, 51030,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 82593, 0, 3,
                                                                       81711, 50310, 81774,
                                                                       29655, 29715, 51120,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 82719, 0, 3,
                                                                       81774, 50355, 81837,
                                                                       29715, 29775, 51210,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 82845, 0, 3,
                                                                       81837, 50400, 81900,
                                                                       29775, 29835, 51300,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 82971, 0, 3,
                                                                       81900, 50445, 81963,
                                                                       29835, 29895, 51390,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 83097, 0, 3,
                                                                       81963, 50490, 82026,
                                                                       29895, 29955, 51480,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 83223, 0, 3,
                                                                       82026, 50535, 82089,
                                                                       29955, 30015, 51570,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 83349, 0, 3,
                                                                       82089, 50580, 82152,
                                                                       30015, 30075, 51660,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 83475, 0, 3,
                                                                       82215, 50850, 82341,
                                                                       30195, 30295, 52050,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 83685, 0, 3,
                                                                       82341, 50940, 82467,
                                                                       30295, 30395, 52200,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 83895, 0, 3,
                                                                       82467, 51030, 82593,
                                                                       30395, 30495, 52350,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 84105, 0, 3,
                                                                       82593, 51120, 82719,
                                                                       30495, 30595, 52500,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 84315, 0, 3,
                                                                       82719, 51210, 82845,
                                                                       30595, 30695, 52650,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 84525, 0, 3,
                                                                       82845, 51300, 82971,
                                                                       30695, 30795, 52800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 84735, 0, 3,
                                                                       82971, 51390, 83097,
                                                                       30795, 30895, 52950,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 84945, 0, 3,
                                                                       83097, 51480, 83223,
                                                                       30895, 30995, 53100,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 85155, 0, 3,
                                                                       83223, 51570, 83349,
                                                                       30995, 31095, 53250,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 85365, 0, 3,
                                                                       83475, 52050, 83685,
                                                                       31295, 31445, 53850,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 85680, 0, 3,
                                                                       83685, 52200, 83895,
                                                                       31445, 31595, 54075,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 85995, 0, 3,
                                                                       83895, 52350, 84105,
                                                                       31595, 31745, 54300,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 86310, 0, 3,
                                                                       84105, 52500, 84315,
                                                                       31745, 31895, 54525,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 86625, 0, 3,
                                                                       84315, 52650, 84525,
                                                                       31895, 32045, 54750,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 86940, 0, 3,
                                                                       84525, 52800, 84735,
                                                                       32045, 32195, 54975,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 87255, 0, 3,
                                                                       84735, 52950, 84945,
                                                                       32195, 32345, 55200,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 87570, 0, 3,
                                                                       84945, 53100, 85155,
                                                                       32345, 32495, 55425,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 87885, 0, 3,
                                                                       85365, 53850, 85680,
                                                                       32795, 33005, 56280,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 88326, 0, 3,
                                                                       85680, 54075, 85995,
                                                                       33005, 33215, 56595,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 88767, 0, 3,
                                                                       85995, 54300, 86310,
                                                                       33215, 33425, 56910,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 89208, 0, 3,
                                                                       86310, 54525, 86625,
                                                                       33425, 33635, 57225,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 89649, 0, 3,
                                                                       86625, 54750, 86940,
                                                                       33635, 33845, 57540,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 90090, 0, 3,
                                                                       86940, 54975, 87255,
                                                                       33845, 34055, 57855,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 90531, 0, 3,
                                                                       87255, 55200, 87570,
                                                                       34055, 34265, 58170,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 90972, 0, 3,
                                                                       87885, 56280, 88326,
                                                                       34685, 34965, 59325,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 91560, 0, 3,
                                                                       88326, 56595, 88767,
                                                                       34965, 35245, 59745,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 92148, 0, 3,
                                                                       88767, 56910, 89208,
                                                                       35245, 35525, 60165,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 92736, 0, 3,
                                                                       89208, 57225, 89649,
                                                                       35525, 35805, 60585,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 93324, 0, 3,
                                                                       89649, 57540, 90090,
                                                                       35805, 36085, 61005,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 93912, 0, 3,
                                                                       90090, 57855, 90531,
                                                                       36085, 36365, 61425,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 94500, 0, 3,
                                                                       90972, 59325, 91560,
                                                                       36925, 37285, 62925,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 95256, 0, 3,
                                                                       91560, 59745, 92148,
                                                                       37285, 37645, 63465,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 96012, 0, 3,
                                                                       92148, 60165, 92736,
                                                                       37645, 38005, 64005,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 96768, 0, 3,
                                                                       92736, 60585, 93324,
                                                                       38005, 38365, 64545,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 97524, 0, 3,
                                                                       93324, 61005, 93912,
                                                                       38365, 38725, 65085,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 98280, 0, 3,
                                                                       94500, 62925, 95256,
                                                                       39445, 39895, 66975,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 99225, 0, 3,
                                                                       95256, 63465, 96012,
                                                                       39895, 40345, 67650,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 100170, 0, 3,
                                                                       96012, 64005, 96768,
                                                                       40345, 40795, 68325,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 101115, 0, 3,
                                                                       96768, 64545, 97524,
                                                                       40795, 41245, 69000,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 102060, 0, 3,
                                                                       98280, 66975, 99225,
                                                                       42145, 42695, 71325,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 103215, 0, 3,
                                                                       99225, 67650, 100170,
                                                                       42695, 43245, 72150,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 104370, 0, 3,
                                                                       100170, 68325, 101115,
                                                                       43245, 43795, 72975,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 105525, 0, 3,
                                                                       102060, 71325, 103215,
                                                                       44895, 45555, 75780,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 106911, 0, 3,
                                                                       103215, 72150, 104370,
                                                                       45555, 46215, 76770,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 108297, 0, 3,
                                                                       105525, 75780, 106911,
                                                                       47535, 48315, 80100,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 109935, 3, 49875,
                                                                       49890, 81270, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 109963, 3, 49890,
                                                                       49905, 81291, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 109991, 3, 49905,
                                                                       49920, 81312, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110019, 3, 49920,
                                                                       49935, 81333, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110047, 3, 49935,
                                                                       49950, 81354, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110075, 3, 49950,
                                                                       49965, 81375, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110103, 3, 49965,
                                                                       49980, 81396, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110131, 3, 49980,
                                                                       49995, 81417, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110159, 3, 49995,
                                                                       50010, 81438, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110187, 3, 50010,
                                                                       50025, 81459, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110215, 3, 50025,
                                                                       50040, 81480, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 110243, 3, 50040,
                                                                       50055, 81501, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 110271, 0, 3,
                                                                       109935, 81270, 109963,
                                                                       50085, 50130, 81522,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 110355, 0, 3,
                                                                       109963, 81291, 109991,
                                                                       50130, 50175, 81585,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 110439, 0, 3,
                                                                       109991, 81312, 110019,
                                                                       50175, 50220, 81648,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 110523, 0, 3,
                                                                       110019, 81333, 110047,
                                                                       50220, 50265, 81711,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 110607, 0, 3,
                                                                       110047, 81354, 110075,
                                                                       50265, 50310, 81774,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 110691, 0, 3,
                                                                       110075, 81375, 110103,
                                                                       50310, 50355, 81837,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 110775, 0, 3,
                                                                       110103, 81396, 110131,
                                                                       50355, 50400, 81900,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 110859, 0, 3,
                                                                       110131, 81417, 110159,
                                                                       50400, 50445, 81963,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 110943, 0, 3,
                                                                       110159, 81438, 110187,
                                                                       50445, 50490, 82026,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 111027, 0, 3,
                                                                       110187, 81459, 110215,
                                                                       50490, 50535, 82089,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 111111, 0, 3,
                                                                       110215, 81480, 110243,
                                                                       50535, 50580, 82152,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111195, 0, 3,
                                                                       110271, 81522, 110355,
                                                                       50670, 50760, 82215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111363, 0, 3,
                                                                       110355, 81585, 110439,
                                                                       50760, 50850, 82341,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111531, 0, 3,
                                                                       110439, 81648, 110523,
                                                                       50850, 50940, 82467,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111699, 0, 3,
                                                                       110523, 81711, 110607,
                                                                       50940, 51030, 82593,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 111867, 0, 3,
                                                                       110607, 81774, 110691,
                                                                       51030, 51120, 82719,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 112035, 0, 3,
                                                                       110691, 81837, 110775,
                                                                       51120, 51210, 82845,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 112203, 0, 3,
                                                                       110775, 81900, 110859,
                                                                       51210, 51300, 82971,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 112371, 0, 3,
                                                                       110859, 81963, 110943,
                                                                       51300, 51390, 83097,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 112539, 0, 3,
                                                                       110943, 82026, 111027,
                                                                       51390, 51480, 83223,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 112707, 0, 3,
                                                                       111027, 82089, 111111,
                                                                       51480, 51570, 83349,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 112875, 0, 3,
                                                                       111195, 82215, 111363,
                                                                       51750, 51900, 83475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 113155, 0, 3,
                                                                       111363, 82341, 111531,
                                                                       51900, 52050, 83685,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 113435, 0, 3,
                                                                       111531, 82467, 111699,
                                                                       52050, 52200, 83895,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 113715, 0, 3,
                                                                       111699, 82593, 111867,
                                                                       52200, 52350, 84105,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 113995, 0, 3,
                                                                       111867, 82719, 112035,
                                                                       52350, 52500, 84315,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 114275, 0, 3,
                                                                       112035, 82845, 112203,
                                                                       52500, 52650, 84525,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 114555, 0, 3,
                                                                       112203, 82971, 112371,
                                                                       52650, 52800, 84735,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 114835, 0, 3,
                                                                       112371, 83097, 112539,
                                                                       52800, 52950, 84945,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 115115, 0, 3,
                                                                       112539, 83223, 112707,
                                                                       52950, 53100, 85155,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 115395, 0, 3,
                                                                       112875, 83475, 113155,
                                                                       53400, 53625, 85365,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 115815, 0, 3,
                                                                       113155, 83685, 113435,
                                                                       53625, 53850, 85680,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 116235, 0, 3,
                                                                       113435, 83895, 113715,
                                                                       53850, 54075, 85995,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 116655, 0, 3,
                                                                       113715, 84105, 113995,
                                                                       54075, 54300, 86310,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 117075, 0, 3,
                                                                       113995, 84315, 114275,
                                                                       54300, 54525, 86625,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 117495, 0, 3,
                                                                       114275, 84525, 114555,
                                                                       54525, 54750, 86940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 117915, 0, 3,
                                                                       114555, 84735, 114835,
                                                                       54750, 54975, 87255,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 118335, 0, 3,
                                                                       114835, 84945, 115115,
                                                                       54975, 55200, 87570,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 118755, 0, 3,
                                                                       115395, 85365, 115815,
                                                                       55650, 55965, 87885,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 119343, 0, 3,
                                                                       115815, 85680, 116235,
                                                                       55965, 56280, 88326,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 119931, 0, 3,
                                                                       116235, 85995, 116655,
                                                                       56280, 56595, 88767,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 120519, 0, 3,
                                                                       116655, 86310, 117075,
                                                                       56595, 56910, 89208,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 121107, 0, 3,
                                                                       117075, 86625, 117495,
                                                                       56910, 57225, 89649,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 121695, 0, 3,
                                                                       117495, 86940, 117915,
                                                                       57225, 57540, 90090,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 122283, 0, 3,
                                                                       117915, 87255, 118335,
                                                                       57540, 57855, 90531,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 122871, 0, 3,
                                                                       118755, 87885, 119343,
                                                                       58485, 58905, 90972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 123655, 0, 3,
                                                                       119343, 88326, 119931,
                                                                       58905, 59325, 91560,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 124439, 0, 3,
                                                                       119931, 88767, 120519,
                                                                       59325, 59745, 92148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 125223, 0, 3,
                                                                       120519, 89208, 121107,
                                                                       59745, 60165, 92736,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 126007, 0, 3,
                                                                       121107, 89649, 121695,
                                                                       60165, 60585, 93324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 126791, 0, 3,
                                                                       121695, 90090, 122283,
                                                                       60585, 61005, 93912,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 127575, 0, 3,
                                                                       122871, 90972, 123655,
                                                                       61845, 62385, 94500,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 128583, 0, 3,
                                                                       123655, 91560, 124439,
                                                                       62385, 62925, 95256,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 129591, 0, 3,
                                                                       124439, 92148, 125223,
                                                                       62925, 63465, 96012,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 130599, 0, 3,
                                                                       125223, 92736, 126007,
                                                                       63465, 64005, 96768,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 131607, 0, 3,
                                                                       126007, 93324, 126791,
                                                                       64005, 64545, 97524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 132615, 0, 3,
                                                                       127575, 94500, 128583,
                                                                       65625, 66300, 98280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 133875, 0, 3,
                                                                       128583, 95256, 129591,
                                                                       66300, 66975, 99225,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 135135, 0, 3,
                                                                       129591, 96012, 130599,
                                                                       66975, 67650, 100170,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 136395, 0, 3,
                                                                       130599, 96768, 131607,
                                                                       67650, 68325, 101115,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 137655, 0, 3,
                                                                       132615, 98280, 133875,
                                                                       69675, 70500, 102060,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 139195, 0, 3,
                                                                       133875, 99225, 135135,
                                                                       70500, 71325, 103215,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 140735, 0, 3,
                                                                       135135, 100170, 136395,
                                                                       71325, 72150, 104370,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 142275, 0, 3,
                                                                       137655, 102060, 139195,
                                                                       73800, 74790, 105525,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 144123, 0, 3,
                                                                       139195, 103215, 140735,
                                                                       74790, 75780, 106911,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 145971, 0, 3,
                                                                       142275, 105525, 144123,
                                                                       77760, 78930, 108297,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 148155, 122871, 784, ncols);

                    simdfunc::contract_primitives(buffer, 149303, 127575, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 150779, 132615, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 152624, 137655, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 154879, 142275, 1848, ncols);

                    simdfunc::contract_primitives(buffer, 157585, 145971, 2184, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 148939, 148155, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 150311, 149303, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 152039, 150779, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 154164, 152624, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 156727, 154879, 66, 1, nmax);

        simdtrf::transform_i_inner(buffer, 159769, 157585, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 160783, 148939, 150311, 13, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 161875, 150311, 152039, 13, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 163279, 152039, 154164, 13, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 165034, 154164, 156727, 13, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 167179, 156727, 159769, 13, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 169753, 160783, 161875, 13, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 171937, 161875, 163279, 13, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 174745, 163279, 165034, 13, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 178255, 165034, 167179, 13, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 182545, 169753, 171937, 13, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 186185, 171937, 174745, 13, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 190865, 174745, 178255, 13, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 196715, 182545, 186185, 13, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 202175, 186185, 190865, 13, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 209195, 196715, 202175, 13, nmax);

        simdtrf::transform_i_inner(buffer, 216839, 209195, 21, 13, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 216839, 169, nmax);
    }

    for (size_t m = 0; m < 1859; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
