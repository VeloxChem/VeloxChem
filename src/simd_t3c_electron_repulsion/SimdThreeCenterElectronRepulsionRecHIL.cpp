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


#include "SimdThreeCenterElectronRepulsionRecHIL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSII.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSML.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSL.hpp"
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
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hil_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hil_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 416914, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2431 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 416914, 319873, 17770, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 19,
                                                             ncols, fj, 6, fq);

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

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 76, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 79, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 82, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 85, 0, 3, 8, 9,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 91, 0, 3, 9, 10,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 97, 0, 3, 10, 11,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 103, 0, 3, 11, 12,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 109, 0, 3, 12, 13,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 115, 0, 3, 13, 14,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 121, 0, 3, 14, 15,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 127, 0, 3, 15, 16,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 133, 0, 3, 16, 17,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 139, 0, 3, 17, 18,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 145, 0, 3, 18, 19,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 151, 0, 3, 19, 20,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 157, 0, 3, 20, 21,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 163, 0, 3, 21, 22,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 169, 0, 3, 22, 23,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 175, 0, 3, 23, 24,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 181, 0, 3, 24, 25,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 187, 0, 3, 25, 26,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 28, 31,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 31, 34,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 34, 37,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 37, 40,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 233, 0, 3, 40, 43,
                                                                       109, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 43, 46,
                                                                       115, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 46, 49,
                                                                       121, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 263, 0, 3, 49, 52,
                                                                       127, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 52, 55,
                                                                       133, 139, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 283, 0, 3, 55, 58,
                                                                       139, 145, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 293, 0, 3, 58, 61,
                                                                       145, 151, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 303, 0, 3, 61, 64,
                                                                       151, 157, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 313, 0, 3, 64, 67,
                                                                       157, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 323, 0, 3, 67, 70,
                                                                       163, 169, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 333, 0, 3, 70, 73,
                                                                       169, 175, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 343, 0, 3, 73, 76,
                                                                       175, 181, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 353, 0, 3, 76, 79,
                                                                       181, 187, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 363, 0, 3, 85, 91,
                                                                       193, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 91, 97,
                                                                       203, 213, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 393, 0, 3, 97,
                                                                       103, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 103,
                                                                       109, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 423, 0, 3, 109,
                                                                       115, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 115,
                                                                       121, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 453, 0, 3, 121,
                                                                       127, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 127,
                                                                       133, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 483, 0, 3, 133,
                                                                       139, 273, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 139,
                                                                       145, 283, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 513, 0, 3, 145,
                                                                       151, 293, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 528, 0, 3, 151,
                                                                       157, 303, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 543, 0, 3, 157,
                                                                       163, 313, 323, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 558, 0, 3, 163,
                                                                       169, 323, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 573, 0, 3, 169,
                                                                       175, 333, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 588, 0, 3, 175,
                                                                       181, 343, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 603, 0, 3, 193,
                                                                       203, 363, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 624, 0, 3, 203,
                                                                       213, 378, 393, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 645, 0, 3, 213,
                                                                       223, 393, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 666, 0, 3, 223,
                                                                       233, 408, 423, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 687, 0, 3, 233,
                                                                       243, 423, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 708, 0, 3, 243,
                                                                       253, 438, 453, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 729, 0, 3, 253,
                                                                       263, 453, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 750, 0, 3, 263,
                                                                       273, 468, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 771, 0, 3, 273,
                                                                       283, 483, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 792, 0, 3, 283,
                                                                       293, 498, 513, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 813, 0, 3, 293,
                                                                       303, 513, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 834, 0, 3, 303,
                                                                       313, 528, 543, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 855, 0, 3, 313,
                                                                       323, 543, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 876, 0, 3, 323,
                                                                       333, 558, 573, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 897, 0, 3, 333,
                                                                       343, 573, 588, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 918, 0, 3, 363,
                                                                       378, 603, 624, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 946, 0, 3, 378,
                                                                       393, 624, 645, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 974, 0, 3, 393,
                                                                       408, 645, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 408,
                                                                       423, 666, 687, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 423,
                                                                       438, 687, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 438,
                                                                       453, 708, 729, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 453,
                                                                       468, 729, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 468,
                                                                       483, 750, 771, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 483,
                                                                       498, 771, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1170, 0, 3, 498,
                                                                       513, 792, 813, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1198, 0, 3, 513,
                                                                       528, 813, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 528,
                                                                       543, 834, 855, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1254, 0, 3, 543,
                                                                       558, 855, 876, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1282, 0, 3, 558,
                                                                       573, 876, 897, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 603,
                                                                       624, 918, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1346, 0, 3, 624,
                                                                       645, 946, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1382, 0, 3, 645,
                                                                       666, 974, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1418, 0, 3, 666,
                                                                       687, 1002, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1454, 0, 3, 687,
                                                                       708, 1030, 1058, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1490, 0, 3, 708,
                                                                       729, 1058, 1086, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1526, 0, 3, 729,
                                                                       750, 1086, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1562, 0, 3, 750,
                                                                       771, 1114, 1142, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1598, 0, 3, 771,
                                                                       792, 1142, 1170, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1634, 0, 3, 792,
                                                                       813, 1170, 1198, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1670, 0, 3, 813,
                                                                       834, 1198, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1706, 0, 3, 834,
                                                                       855, 1226, 1254, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1742, 0, 3, 855,
                                                                       876, 1254, 1282, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1778, 0, 3, 918,
                                                                       946, 1310, 1346, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1823, 0, 3, 946,
                                                                       974, 1346, 1382, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1868, 0, 3, 974,
                                                                       1002, 1382, 1418, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1913, 0, 3, 1002,
                                                                       1030, 1418, 1454, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1958, 0, 3, 1030,
                                                                       1058, 1454, 1490, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2003, 0, 3, 1058,
                                                                       1086, 1490, 1526, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2048, 0, 3, 1086,
                                                                       1114, 1526, 1562, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2093, 0, 3, 1114,
                                                                       1142, 1562, 1598, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2138, 0, 3, 1142,
                                                                       1170, 1598, 1634, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2183, 0, 3, 1170,
                                                                       1198, 1634, 1670, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2228, 0, 3, 1198,
                                                                       1226, 1670, 1706, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 2273, 0, 3, 1226,
                                                                       1254, 1706, 1742, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2318, 0, 3, 1310,
                                                                       1346, 1778, 1823, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2373, 0, 3, 1346,
                                                                       1382, 1823, 1868, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2428, 0, 3, 1382,
                                                                       1418, 1868, 1913, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2483, 0, 3, 1418,
                                                                       1454, 1913, 1958, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2538, 0, 3, 1454,
                                                                       1490, 1958, 2003, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2593, 0, 3, 1490,
                                                                       1526, 2003, 2048, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 1526,
                                                                       1562, 2048, 2093, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2703, 0, 3, 1562,
                                                                       1598, 2093, 2138, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2758, 0, 3, 1598,
                                                                       1634, 2138, 2183, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2813, 0, 3, 1634,
                                                                       1670, 2183, 2228, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2868, 0, 3, 1670,
                                                                       1706, 2228, 2273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2923, 0, 3, 1778,
                                                                       1823, 2318, 2373, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2989, 0, 3, 1823,
                                                                       1868, 2373, 2428, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3055, 0, 3, 1868,
                                                                       1913, 2428, 2483, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3121, 0, 3, 1913,
                                                                       1958, 2483, 2538, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3187, 0, 3, 1958,
                                                                       2003, 2538, 2593, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3253, 0, 3, 2003,
                                                                       2048, 2593, 2648, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3319, 0, 3, 2048,
                                                                       2093, 2648, 2703, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3385, 0, 3, 2093,
                                                                       2138, 2703, 2758, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3451, 0, 3, 2138,
                                                                       2183, 2758, 2813, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 3517, 0, 3, 2183,
                                                                       2228, 2813, 2868, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3583, 0, 3, 2318,
                                                                       2373, 2923, 2989, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3661, 0, 3, 2373,
                                                                       2428, 2989, 3055, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3739, 0, 3, 2428,
                                                                       2483, 3055, 3121, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3817, 0, 3, 2483,
                                                                       2538, 3121, 3187, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3895, 0, 3, 2538,
                                                                       2593, 3187, 3253, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 3973, 0, 3, 2593,
                                                                       2648, 3253, 3319, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 4051, 0, 3, 2648,
                                                                       2703, 3319, 3385, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 4129, 0, 3, 2703,
                                                                       2758, 3385, 3451, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 4207, 0, 3, 2758,
                                                                       2813, 3451, 3517, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4285, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4288, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4291, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4294, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4297, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4300, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4303, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4306, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4309, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4312, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4315, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4318, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4321, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4324, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4327, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4330, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4333, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 4336, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4339, 3, 10, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4348, 3, 11, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4357, 3, 12, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4366, 3, 13, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4375, 3, 14, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4384, 3, 15, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4393, 3, 16, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4402, 3, 17, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4411, 3, 18, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4420, 3, 19, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4429, 3, 20, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4438, 3, 21, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4447, 3, 22, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4456, 3, 23, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4465, 3, 24, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4474, 3, 25, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 4483, 3, 26, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4492, 3, 34, 97,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4510, 3, 37, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4528, 3, 40, 109,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4546, 3, 43, 115,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4564, 3, 46, 121,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4582, 3, 49, 127,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4600, 3, 52, 133,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4618, 3, 55, 139,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4636, 3, 58, 145,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4654, 3, 61, 151,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4672, 3, 64, 157,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4690, 3, 67, 163,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4708, 3, 70, 169,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4726, 3, 73, 175,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4744, 3, 76, 181,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 4762, 3, 79, 187,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4780, 3, 97, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4810, 3, 103, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4840, 3, 109, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4870, 3, 115, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4900, 3, 121, 253,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4930, 3, 127, 263,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4960, 3, 133, 273,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 4990, 3, 139, 283,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5020, 3, 145, 293,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5050, 3, 151, 303,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5080, 3, 157, 313,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5110, 3, 163, 323,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5140, 3, 169, 333,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5170, 3, 175, 343,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 5200, 3, 181, 353,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5230, 3, 213, 393,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5275, 3, 223, 408,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5320, 3, 233, 423,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5365, 3, 243, 438,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5410, 3, 253, 453,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5455, 3, 263, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5500, 3, 273, 483,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5545, 3, 283, 498,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5590, 3, 293, 513,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5635, 3, 303, 528,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5680, 3, 313, 543,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5725, 3, 323, 558,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5770, 3, 333, 573,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 5815, 3, 343, 588,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5860, 3, 393, 645,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5923, 3, 408, 666,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 5986, 3, 423, 687,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6049, 3, 438, 708,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6112, 3, 453, 729,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6175, 3, 468, 750,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6238, 3, 483, 771,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6301, 3, 498, 792,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6364, 3, 513, 813,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6427, 3, 528, 834,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6490, 3, 543, 855,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6553, 3, 558, 876,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 6616, 3, 573, 897,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6679, 3, 645, 974,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6763, 3, 666,
                                                                       1002, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6847, 3, 687,
                                                                       1030, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 6931, 3, 708,
                                                                       1058, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7015, 3, 729,
                                                                       1086, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7099, 3, 750,
                                                                       1114, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7183, 3, 771,
                                                                       1142, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7267, 3, 792,
                                                                       1170, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7351, 3, 813,
                                                                       1198, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7435, 3, 834,
                                                                       1226, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7519, 3, 855,
                                                                       1254, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 7603, 3, 876,
                                                                       1282, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7687, 3, 974,
                                                                       1382, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7795, 3, 1002,
                                                                       1418, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 7903, 3, 1030,
                                                                       1454, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8011, 3, 1058,
                                                                       1490, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8119, 3, 1086,
                                                                       1526, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8227, 3, 1114,
                                                                       1562, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8335, 3, 1142,
                                                                       1598, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8443, 3, 1170,
                                                                       1634, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8551, 3, 1198,
                                                                       1670, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8659, 3, 1226,
                                                                       1706, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 8767, 3, 1254,
                                                                       1742, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8875, 3, 1382,
                                                                       1868, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9010, 3, 1418,
                                                                       1913, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9145, 3, 1454,
                                                                       1958, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9280, 3, 1490,
                                                                       2003, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9415, 3, 1526,
                                                                       2048, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9550, 3, 1562,
                                                                       2093, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9685, 3, 1598,
                                                                       2138, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9820, 3, 1634,
                                                                       2183, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 9955, 3, 1670,
                                                                       2228, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 10090, 3, 1706,
                                                                       2273, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10225, 3, 1868,
                                                                       2428, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10390, 3, 1913,
                                                                       2483, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10555, 3, 1958,
                                                                       2538, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10720, 3, 2003,
                                                                       2593, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 10885, 3, 2048,
                                                                       2648, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11050, 3, 2093,
                                                                       2703, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11215, 3, 2138,
                                                                       2758, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11380, 3, 2183,
                                                                       2813, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 11545, 3, 2228,
                                                                       2868, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11710, 3, 2428,
                                                                       3055, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 11908, 3, 2483,
                                                                       3121, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12106, 3, 2538,
                                                                       3187, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12304, 3, 2593,
                                                                       3253, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12502, 3, 2648,
                                                                       3319, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12700, 3, 2703,
                                                                       3385, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 12898, 3, 2758,
                                                                       3451, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 13096, 3, 2813,
                                                                       3517, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13294, 3, 3055,
                                                                       3739, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13528, 3, 3121,
                                                                       3817, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13762, 3, 3187,
                                                                       3895, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 13996, 3, 3253,
                                                                       3973, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 14230, 3, 3319,
                                                                       4051, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 14464, 3, 3385,
                                                                       4129, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 14698, 3, 3451,
                                                                       4207, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14932, 3, 8, 9,
                                                                       4285, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14938, 3, 9, 10,
                                                                       4288, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14944, 3, 10, 11,
                                                                       4291, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14950, 3, 11, 12,
                                                                       4294, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14956, 3, 12, 13,
                                                                       4297, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14962, 3, 13, 14,
                                                                       4300, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14968, 3, 14, 15,
                                                                       4303, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14974, 3, 15, 16,
                                                                       4306, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14980, 3, 16, 17,
                                                                       4309, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14986, 3, 17, 18,
                                                                       4312, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14992, 3, 18, 19,
                                                                       4315, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 14998, 3, 19, 20,
                                                                       4318, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15004, 3, 20, 21,
                                                                       4321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15010, 3, 21, 22,
                                                                       4324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15016, 3, 22, 23,
                                                                       4327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15022, 3, 23, 24,
                                                                       4330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15028, 3, 24, 25,
                                                                       4333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 15034, 3, 25, 26,
                                                                       4336, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15040, 0, 3,
                                                                       14932, 4285, 14938, 4339,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15058, 0, 3,
                                                                       14938, 4288, 14944, 4348,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15076, 0, 3,
                                                                       14944, 4291, 14950, 4357,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15094, 0, 3,
                                                                       14950, 4294, 14956, 4366,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15112, 0, 3,
                                                                       14956, 4297, 14962, 4375,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15130, 0, 3,
                                                                       14962, 4300, 14968, 4384,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15148, 0, 3,
                                                                       14968, 4303, 14974, 4393,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15166, 0, 3,
                                                                       14974, 4306, 14980, 4402,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15184, 0, 3,
                                                                       14980, 4309, 14986, 4411,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15202, 0, 3,
                                                                       14986, 4312, 14992, 4420,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15220, 0, 3,
                                                                       14992, 4315, 14998, 4429,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15238, 0, 3,
                                                                       14998, 4318, 15004, 4438,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15256, 0, 3,
                                                                       15004, 4321, 15010, 4447,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15274, 0, 3,
                                                                       15010, 4324, 15016, 4456,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15292, 0, 3,
                                                                       15016, 4327, 15022, 4465,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15310, 0, 3,
                                                                       15022, 4330, 15028, 4474,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 15328, 0, 3,
                                                                       15028, 4333, 15034, 4483,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15346, 0, 3,
                                                                       15040, 4339, 15058, 85,
                                                                       91, 4492, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15382, 0, 3,
                                                                       15058, 4348, 15076, 91,
                                                                       97, 4510, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15418, 0, 3,
                                                                       15076, 4357, 15094, 97,
                                                                       103, 4528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15454, 0, 3,
                                                                       15094, 4366, 15112, 103,
                                                                       109, 4546, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15490, 0, 3,
                                                                       15112, 4375, 15130, 109,
                                                                       115, 4564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15526, 0, 3,
                                                                       15130, 4384, 15148, 115,
                                                                       121, 4582, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15562, 0, 3,
                                                                       15148, 4393, 15166, 121,
                                                                       127, 4600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15598, 0, 3,
                                                                       15166, 4402, 15184, 127,
                                                                       133, 4618, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15634, 0, 3,
                                                                       15184, 4411, 15202, 133,
                                                                       139, 4636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15670, 0, 3,
                                                                       15202, 4420, 15220, 139,
                                                                       145, 4654, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15706, 0, 3,
                                                                       15220, 4429, 15238, 145,
                                                                       151, 4672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15742, 0, 3,
                                                                       15238, 4438, 15256, 151,
                                                                       157, 4690, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15778, 0, 3,
                                                                       15256, 4447, 15274, 157,
                                                                       163, 4708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15814, 0, 3,
                                                                       15274, 4456, 15292, 163,
                                                                       169, 4726, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15850, 0, 3,
                                                                       15292, 4465, 15310, 169,
                                                                       175, 4744, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 15886, 0, 3,
                                                                       15310, 4474, 15328, 175,
                                                                       181, 4762, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15922, 0, 3,
                                                                       15346, 4492, 15382, 193,
                                                                       203, 4780, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 15982, 0, 3,
                                                                       15382, 4510, 15418, 203,
                                                                       213, 4810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16042, 0, 3,
                                                                       15418, 4528, 15454, 213,
                                                                       223, 4840, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16102, 0, 3,
                                                                       15454, 4546, 15490, 223,
                                                                       233, 4870, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16162, 0, 3,
                                                                       15490, 4564, 15526, 233,
                                                                       243, 4900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16222, 0, 3,
                                                                       15526, 4582, 15562, 243,
                                                                       253, 4930, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16282, 0, 3,
                                                                       15562, 4600, 15598, 253,
                                                                       263, 4960, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16342, 0, 3,
                                                                       15598, 4618, 15634, 263,
                                                                       273, 4990, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16402, 0, 3,
                                                                       15634, 4636, 15670, 273,
                                                                       283, 5020, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16462, 0, 3,
                                                                       15670, 4654, 15706, 283,
                                                                       293, 5050, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16522, 0, 3,
                                                                       15706, 4672, 15742, 293,
                                                                       303, 5080, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16582, 0, 3,
                                                                       15742, 4690, 15778, 303,
                                                                       313, 5110, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16642, 0, 3,
                                                                       15778, 4708, 15814, 313,
                                                                       323, 5140, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16702, 0, 3,
                                                                       15814, 4726, 15850, 323,
                                                                       333, 5170, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 16762, 0, 3,
                                                                       15850, 4744, 15886, 333,
                                                                       343, 5200, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16822, 0, 3,
                                                                       15922, 4780, 15982, 363,
                                                                       378, 5230, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 16912, 0, 3,
                                                                       15982, 4810, 16042, 378,
                                                                       393, 5275, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17002, 0, 3,
                                                                       16042, 4840, 16102, 393,
                                                                       408, 5320, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17092, 0, 3,
                                                                       16102, 4870, 16162, 408,
                                                                       423, 5365, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17182, 0, 3,
                                                                       16162, 4900, 16222, 423,
                                                                       438, 5410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17272, 0, 3,
                                                                       16222, 4930, 16282, 438,
                                                                       453, 5455, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17362, 0, 3,
                                                                       16282, 4960, 16342, 453,
                                                                       468, 5500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17452, 0, 3,
                                                                       16342, 4990, 16402, 468,
                                                                       483, 5545, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17542, 0, 3,
                                                                       16402, 5020, 16462, 483,
                                                                       498, 5590, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17632, 0, 3,
                                                                       16462, 5050, 16522, 498,
                                                                       513, 5635, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17722, 0, 3,
                                                                       16522, 5080, 16582, 513,
                                                                       528, 5680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17812, 0, 3,
                                                                       16582, 5110, 16642, 528,
                                                                       543, 5725, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17902, 0, 3,
                                                                       16642, 5140, 16702, 543,
                                                                       558, 5770, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 17992, 0, 3,
                                                                       16702, 5170, 16762, 558,
                                                                       573, 5815, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18082, 0, 3,
                                                                       16822, 5230, 16912, 603,
                                                                       624, 5860, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18208, 0, 3,
                                                                       16912, 5275, 17002, 624,
                                                                       645, 5923, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18334, 0, 3,
                                                                       17002, 5320, 17092, 645,
                                                                       666, 5986, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18460, 0, 3,
                                                                       17092, 5365, 17182, 666,
                                                                       687, 6049, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18586, 0, 3,
                                                                       17182, 5410, 17272, 687,
                                                                       708, 6112, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18712, 0, 3,
                                                                       17272, 5455, 17362, 708,
                                                                       729, 6175, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18838, 0, 3,
                                                                       17362, 5500, 17452, 729,
                                                                       750, 6238, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 18964, 0, 3,
                                                                       17452, 5545, 17542, 750,
                                                                       771, 6301, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19090, 0, 3,
                                                                       17542, 5590, 17632, 771,
                                                                       792, 6364, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19216, 0, 3,
                                                                       17632, 5635, 17722, 792,
                                                                       813, 6427, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19342, 0, 3,
                                                                       17722, 5680, 17812, 813,
                                                                       834, 6490, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19468, 0, 3,
                                                                       17812, 5725, 17902, 834,
                                                                       855, 6553, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 19594, 0, 3,
                                                                       17902, 5770, 17992, 855,
                                                                       876, 6616, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19720, 0, 3,
                                                                       18082, 5860, 18208, 918,
                                                                       946, 6679, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 19888, 0, 3,
                                                                       18208, 5923, 18334, 946,
                                                                       974, 6763, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20056, 0, 3,
                                                                       18334, 5986, 18460, 974,
                                                                       1002, 6847, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20224, 0, 3,
                                                                       18460, 6049, 18586, 1002,
                                                                       1030, 6931, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20392, 0, 3,
                                                                       18586, 6112, 18712, 1030,
                                                                       1058, 7015, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20560, 0, 3,
                                                                       18712, 6175, 18838, 1058,
                                                                       1086, 7099, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20728, 0, 3,
                                                                       18838, 6238, 18964, 1086,
                                                                       1114, 7183, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 20896, 0, 3,
                                                                       18964, 6301, 19090, 1114,
                                                                       1142, 7267, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21064, 0, 3,
                                                                       19090, 6364, 19216, 1142,
                                                                       1170, 7351, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21232, 0, 3,
                                                                       19216, 6427, 19342, 1170,
                                                                       1198, 7435, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21400, 0, 3,
                                                                       19342, 6490, 19468, 1198,
                                                                       1226, 7519, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 21568, 0, 3,
                                                                       19468, 6553, 19594, 1226,
                                                                       1254, 7603, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21736, 0, 3,
                                                                       19720, 6679, 19888, 1310,
                                                                       1346, 7687, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 21952, 0, 3,
                                                                       19888, 6763, 20056, 1346,
                                                                       1382, 7795, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22168, 0, 3,
                                                                       20056, 6847, 20224, 1382,
                                                                       1418, 7903, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22384, 0, 3,
                                                                       20224, 6931, 20392, 1418,
                                                                       1454, 8011, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22600, 0, 3,
                                                                       20392, 7015, 20560, 1454,
                                                                       1490, 8119, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 22816, 0, 3,
                                                                       20560, 7099, 20728, 1490,
                                                                       1526, 8227, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23032, 0, 3,
                                                                       20728, 7183, 20896, 1526,
                                                                       1562, 8335, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23248, 0, 3,
                                                                       20896, 7267, 21064, 1562,
                                                                       1598, 8443, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23464, 0, 3,
                                                                       21064, 7351, 21232, 1598,
                                                                       1634, 8551, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23680, 0, 3,
                                                                       21232, 7435, 21400, 1634,
                                                                       1670, 8659, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 23896, 0, 3,
                                                                       21400, 7519, 21568, 1670,
                                                                       1706, 8767, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24112, 0, 3,
                                                                       21736, 7687, 21952, 1778,
                                                                       1823, 8875, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24382, 0, 3,
                                                                       21952, 7795, 22168, 1823,
                                                                       1868, 9010, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24652, 0, 3,
                                                                       22168, 7903, 22384, 1868,
                                                                       1913, 9145, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 24922, 0, 3,
                                                                       22384, 8011, 22600, 1913,
                                                                       1958, 9280, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 25192, 0, 3,
                                                                       22600, 8119, 22816, 1958,
                                                                       2003, 9415, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 25462, 0, 3,
                                                                       22816, 8227, 23032, 2003,
                                                                       2048, 9550, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 25732, 0, 3,
                                                                       23032, 8335, 23248, 2048,
                                                                       2093, 9685, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26002, 0, 3,
                                                                       23248, 8443, 23464, 2093,
                                                                       2138, 9820, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26272, 0, 3,
                                                                       23464, 8551, 23680, 2138,
                                                                       2183, 9955, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 26542, 0, 3,
                                                                       23680, 8659, 23896, 2183,
                                                                       2228, 10090, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 26812, 0, 3,
                                                                       24112, 8875, 24382, 2318,
                                                                       2373, 10225, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 27142, 0, 3,
                                                                       24382, 9010, 24652, 2373,
                                                                       2428, 10390, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 27472, 0, 3,
                                                                       24652, 9145, 24922, 2428,
                                                                       2483, 10555, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 27802, 0, 3,
                                                                       24922, 9280, 25192, 2483,
                                                                       2538, 10720, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 28132, 0, 3,
                                                                       25192, 9415, 25462, 2538,
                                                                       2593, 10885, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 28462, 0, 3,
                                                                       25462, 9550, 25732, 2593,
                                                                       2648, 11050, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 28792, 0, 3,
                                                                       25732, 9685, 26002, 2648,
                                                                       2703, 11215, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 29122, 0, 3,
                                                                       26002, 9820, 26272, 2703,
                                                                       2758, 11380, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 29452, 0, 3,
                                                                       26272, 9955, 26542, 2758,
                                                                       2813, 11545, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 29782, 0, 3,
                                                                       26812, 10225, 27142, 2923,
                                                                       2989, 11710, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 30178, 0, 3,
                                                                       27142, 10390, 27472, 2989,
                                                                       3055, 11908, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 30574, 0, 3,
                                                                       27472, 10555, 27802, 3055,
                                                                       3121, 12106, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 30970, 0, 3,
                                                                       27802, 10720, 28132, 3121,
                                                                       3187, 12304, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 31366, 0, 3,
                                                                       28132, 10885, 28462, 3187,
                                                                       3253, 12502, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 31762, 0, 3,
                                                                       28462, 11050, 28792, 3253,
                                                                       3319, 12700, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 32158, 0, 3,
                                                                       28792, 11215, 29122, 3319,
                                                                       3385, 12898, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 32554, 0, 3,
                                                                       29122, 11380, 29452, 3385,
                                                                       3451, 13096, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 32950, 0, 3,
                                                                       29782, 11710, 30178, 3583,
                                                                       3661, 13294, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 33418, 0, 3,
                                                                       30178, 11908, 30574, 3661,
                                                                       3739, 13528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 33886, 0, 3,
                                                                       30574, 12106, 30970, 3739,
                                                                       3817, 13762, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 34354, 0, 3,
                                                                       30970, 12304, 31366, 3817,
                                                                       3895, 13996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 34822, 0, 3,
                                                                       31366, 12502, 31762, 3895,
                                                                       3973, 14230, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 35290, 0, 3,
                                                                       31762, 12700, 32158, 3973,
                                                                       4051, 14464, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 35758, 0, 3,
                                                                       32158, 12898, 32554, 4051,
                                                                       4129, 14698, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36226, 3, 4285,
                                                                       4288, 14944, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36236, 3, 4288,
                                                                       4291, 14950, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36246, 3, 4291,
                                                                       4294, 14956, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36256, 3, 4294,
                                                                       4297, 14962, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36266, 3, 4297,
                                                                       4300, 14968, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36276, 3, 4300,
                                                                       4303, 14974, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36286, 3, 4303,
                                                                       4306, 14980, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36296, 3, 4306,
                                                                       4309, 14986, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36306, 3, 4309,
                                                                       4312, 14992, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36316, 3, 4312,
                                                                       4315, 14998, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36326, 3, 4315,
                                                                       4318, 15004, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36336, 3, 4318,
                                                                       4321, 15010, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36346, 3, 4321,
                                                                       4324, 15016, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36356, 3, 4324,
                                                                       4327, 15022, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36366, 3, 4327,
                                                                       4330, 15028, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 36376, 3, 4330,
                                                                       4333, 15034, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36386, 0, 3,
                                                                       36226, 14944, 36236,
                                                                       15076, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36416, 0, 3,
                                                                       36236, 14950, 36246,
                                                                       15094, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36446, 0, 3,
                                                                       36246, 14956, 36256,
                                                                       15112, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36476, 0, 3,
                                                                       36256, 14962, 36266,
                                                                       15130, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36506, 0, 3,
                                                                       36266, 14968, 36276,
                                                                       15148, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36536, 0, 3,
                                                                       36276, 14974, 36286,
                                                                       15166, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36566, 0, 3,
                                                                       36286, 14980, 36296,
                                                                       15184, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36596, 0, 3,
                                                                       36296, 14986, 36306,
                                                                       15202, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36626, 0, 3,
                                                                       36306, 14992, 36316,
                                                                       15220, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36656, 0, 3,
                                                                       36316, 14998, 36326,
                                                                       15238, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36686, 0, 3,
                                                                       36326, 15004, 36336,
                                                                       15256, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36716, 0, 3,
                                                                       36336, 15010, 36346,
                                                                       15274, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36746, 0, 3,
                                                                       36346, 15016, 36356,
                                                                       15292, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36776, 0, 3,
                                                                       36356, 15022, 36366,
                                                                       15310, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 36806, 0, 3,
                                                                       36366, 15028, 36376,
                                                                       15328, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36836, 0, 3,
                                                                       36386, 15076, 36416, 4492,
                                                                       4510, 15418, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36896, 0, 3,
                                                                       36416, 15094, 36446, 4510,
                                                                       4528, 15454, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 36956, 0, 3,
                                                                       36446, 15112, 36476, 4528,
                                                                       4546, 15490, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37016, 0, 3,
                                                                       36476, 15130, 36506, 4546,
                                                                       4564, 15526, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37076, 0, 3,
                                                                       36506, 15148, 36536, 4564,
                                                                       4582, 15562, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37136, 0, 3,
                                                                       36536, 15166, 36566, 4582,
                                                                       4600, 15598, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37196, 0, 3,
                                                                       36566, 15184, 36596, 4600,
                                                                       4618, 15634, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37256, 0, 3,
                                                                       36596, 15202, 36626, 4618,
                                                                       4636, 15670, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37316, 0, 3,
                                                                       36626, 15220, 36656, 4636,
                                                                       4654, 15706, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37376, 0, 3,
                                                                       36656, 15238, 36686, 4654,
                                                                       4672, 15742, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37436, 0, 3,
                                                                       36686, 15256, 36716, 4672,
                                                                       4690, 15778, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37496, 0, 3,
                                                                       36716, 15274, 36746, 4690,
                                                                       4708, 15814, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37556, 0, 3,
                                                                       36746, 15292, 36776, 4708,
                                                                       4726, 15850, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 37616, 0, 3,
                                                                       36776, 15310, 36806, 4726,
                                                                       4744, 15886, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 37676, 0, 3,
                                                                       36836, 15418, 36896, 4780,
                                                                       4810, 16042, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 37776, 0, 3,
                                                                       36896, 15454, 36956, 4810,
                                                                       4840, 16102, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 37876, 0, 3,
                                                                       36956, 15490, 37016, 4840,
                                                                       4870, 16162, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 37976, 0, 3,
                                                                       37016, 15526, 37076, 4870,
                                                                       4900, 16222, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38076, 0, 3,
                                                                       37076, 15562, 37136, 4900,
                                                                       4930, 16282, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38176, 0, 3,
                                                                       37136, 15598, 37196, 4930,
                                                                       4960, 16342, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38276, 0, 3,
                                                                       37196, 15634, 37256, 4960,
                                                                       4990, 16402, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38376, 0, 3,
                                                                       37256, 15670, 37316, 4990,
                                                                       5020, 16462, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38476, 0, 3,
                                                                       37316, 15706, 37376, 5020,
                                                                       5050, 16522, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38576, 0, 3,
                                                                       37376, 15742, 37436, 5050,
                                                                       5080, 16582, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38676, 0, 3,
                                                                       37436, 15778, 37496, 5080,
                                                                       5110, 16642, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38776, 0, 3,
                                                                       37496, 15814, 37556, 5110,
                                                                       5140, 16702, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 38876, 0, 3,
                                                                       37556, 15850, 37616, 5140,
                                                                       5170, 16762, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 38976, 0, 3,
                                                                       37676, 16042, 37776, 5230,
                                                                       5275, 17002, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39126, 0, 3,
                                                                       37776, 16102, 37876, 5275,
                                                                       5320, 17092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39276, 0, 3,
                                                                       37876, 16162, 37976, 5320,
                                                                       5365, 17182, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39426, 0, 3,
                                                                       37976, 16222, 38076, 5365,
                                                                       5410, 17272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39576, 0, 3,
                                                                       38076, 16282, 38176, 5410,
                                                                       5455, 17362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39726, 0, 3,
                                                                       38176, 16342, 38276, 5455,
                                                                       5500, 17452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 39876, 0, 3,
                                                                       38276, 16402, 38376, 5500,
                                                                       5545, 17542, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40026, 0, 3,
                                                                       38376, 16462, 38476, 5545,
                                                                       5590, 17632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40176, 0, 3,
                                                                       38476, 16522, 38576, 5590,
                                                                       5635, 17722, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40326, 0, 3,
                                                                       38576, 16582, 38676, 5635,
                                                                       5680, 17812, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40476, 0, 3,
                                                                       38676, 16642, 38776, 5680,
                                                                       5725, 17902, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 40626, 0, 3,
                                                                       38776, 16702, 38876, 5725,
                                                                       5770, 17992, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40776, 0, 3,
                                                                       38976, 17002, 39126, 5860,
                                                                       5923, 18334, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 40986, 0, 3,
                                                                       39126, 17092, 39276, 5923,
                                                                       5986, 18460, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41196, 0, 3,
                                                                       39276, 17182, 39426, 5986,
                                                                       6049, 18586, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41406, 0, 3,
                                                                       39426, 17272, 39576, 6049,
                                                                       6112, 18712, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41616, 0, 3,
                                                                       39576, 17362, 39726, 6112,
                                                                       6175, 18838, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 41826, 0, 3,
                                                                       39726, 17452, 39876, 6175,
                                                                       6238, 18964, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42036, 0, 3,
                                                                       39876, 17542, 40026, 6238,
                                                                       6301, 19090, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42246, 0, 3,
                                                                       40026, 17632, 40176, 6301,
                                                                       6364, 19216, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42456, 0, 3,
                                                                       40176, 17722, 40326, 6364,
                                                                       6427, 19342, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42666, 0, 3,
                                                                       40326, 17812, 40476, 6427,
                                                                       6490, 19468, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 42876, 0, 3,
                                                                       40476, 17902, 40626, 6490,
                                                                       6553, 19594, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43086, 0, 3,
                                                                       40776, 18334, 40986, 6679,
                                                                       6763, 20056, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43366, 0, 3,
                                                                       40986, 18460, 41196, 6763,
                                                                       6847, 20224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43646, 0, 3,
                                                                       41196, 18586, 41406, 6847,
                                                                       6931, 20392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 43926, 0, 3,
                                                                       41406, 18712, 41616, 6931,
                                                                       7015, 20560, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44206, 0, 3,
                                                                       41616, 18838, 41826, 7015,
                                                                       7099, 20728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44486, 0, 3,
                                                                       41826, 18964, 42036, 7099,
                                                                       7183, 20896, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 44766, 0, 3,
                                                                       42036, 19090, 42246, 7183,
                                                                       7267, 21064, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 45046, 0, 3,
                                                                       42246, 19216, 42456, 7267,
                                                                       7351, 21232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 45326, 0, 3,
                                                                       42456, 19342, 42666, 7351,
                                                                       7435, 21400, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 45606, 0, 3,
                                                                       42666, 19468, 42876, 7435,
                                                                       7519, 21568, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 45886, 0, 3,
                                                                       43086, 20056, 43366, 7687,
                                                                       7795, 22168, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46246, 0, 3,
                                                                       43366, 20224, 43646, 7795,
                                                                       7903, 22384, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46606, 0, 3,
                                                                       43646, 20392, 43926, 7903,
                                                                       8011, 22600, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 46966, 0, 3,
                                                                       43926, 20560, 44206, 8011,
                                                                       8119, 22816, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 47326, 0, 3,
                                                                       44206, 20728, 44486, 8119,
                                                                       8227, 23032, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 47686, 0, 3,
                                                                       44486, 20896, 44766, 8227,
                                                                       8335, 23248, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48046, 0, 3,
                                                                       44766, 21064, 45046, 8335,
                                                                       8443, 23464, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48406, 0, 3,
                                                                       45046, 21232, 45326, 8443,
                                                                       8551, 23680, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 48766, 0, 3,
                                                                       45326, 21400, 45606, 8551,
                                                                       8659, 23896, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 49126, 0, 3,
                                                                       45886, 22168, 46246, 8875,
                                                                       9010, 24652, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 49576, 0, 3,
                                                                       46246, 22384, 46606, 9010,
                                                                       9145, 24922, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 50026, 0, 3,
                                                                       46606, 22600, 46966, 9145,
                                                                       9280, 25192, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 50476, 0, 3,
                                                                       46966, 22816, 47326, 9280,
                                                                       9415, 25462, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 50926, 0, 3,
                                                                       47326, 23032, 47686, 9415,
                                                                       9550, 25732, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 51376, 0, 3,
                                                                       47686, 23248, 48046, 9550,
                                                                       9685, 26002, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 51826, 0, 3,
                                                                       48046, 23464, 48406, 9685,
                                                                       9820, 26272, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 52276, 0, 3,
                                                                       48406, 23680, 48766, 9820,
                                                                       9955, 26542, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 52726, 0, 3,
                                                                       49126, 24652, 49576,
                                                                       10225, 10390, 27472,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 53276, 0, 3,
                                                                       49576, 24922, 50026,
                                                                       10390, 10555, 27802,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 53826, 0, 3,
                                                                       50026, 25192, 50476,
                                                                       10555, 10720, 28132,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 54376, 0, 3,
                                                                       50476, 25462, 50926,
                                                                       10720, 10885, 28462,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 54926, 0, 3,
                                                                       50926, 25732, 51376,
                                                                       10885, 11050, 28792,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 55476, 0, 3,
                                                                       51376, 26002, 51826,
                                                                       11050, 11215, 29122,
                                                                       ncols, gamma, p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 56026, 0, 3,
                                                                       51826, 26272, 52276,
                                                                       11215, 11380, 29452,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 56576, 0, 3,
                                                                       52726, 27472, 53276,
                                                                       11710, 11908, 30574,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 57236, 0, 3,
                                                                       53276, 27802, 53826,
                                                                       11908, 12106, 30970,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 57896, 0, 3,
                                                                       53826, 28132, 54376,
                                                                       12106, 12304, 31366,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 58556, 0, 3,
                                                                       54376, 28462, 54926,
                                                                       12304, 12502, 31762,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 59216, 0, 3,
                                                                       54926, 28792, 55476,
                                                                       12502, 12700, 32158,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 59876, 0, 3,
                                                                       55476, 29122, 56026,
                                                                       12700, 12898, 32554,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 60536, 0, 3,
                                                                       56576, 30574, 57236,
                                                                       13294, 13528, 33886,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 61316, 0, 3,
                                                                       57236, 30970, 57896,
                                                                       13528, 13762, 34354,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 62096, 0, 3,
                                                                       57896, 31366, 58556,
                                                                       13762, 13996, 34822,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 62876, 0, 3,
                                                                       58556, 31762, 59216,
                                                                       13996, 14230, 35290,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 63656, 0, 3,
                                                                       59216, 32158, 59876,
                                                                       14230, 14464, 35758,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64436, 3, 14932,
                                                                       14938, 36226, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64451, 3, 14938,
                                                                       14944, 36236, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64466, 3, 14944,
                                                                       14950, 36246, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64481, 3, 14950,
                                                                       14956, 36256, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64496, 3, 14956,
                                                                       14962, 36266, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64511, 3, 14962,
                                                                       14968, 36276, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64526, 3, 14968,
                                                                       14974, 36286, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64541, 3, 14974,
                                                                       14980, 36296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64556, 3, 14980,
                                                                       14986, 36306, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64571, 3, 14986,
                                                                       14992, 36316, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64586, 3, 14992,
                                                                       14998, 36326, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64601, 3, 14998,
                                                                       15004, 36336, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64616, 3, 15004,
                                                                       15010, 36346, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64631, 3, 15010,
                                                                       15016, 36356, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64646, 3, 15016,
                                                                       15022, 36366, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 64661, 3, 15022,
                                                                       15028, 36376, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64676, 0, 3,
                                                                       64436, 36226, 64451,
                                                                       15040, 15058, 36386,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64721, 0, 3,
                                                                       64451, 36236, 64466,
                                                                       15058, 15076, 36416,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64766, 0, 3,
                                                                       64466, 36246, 64481,
                                                                       15076, 15094, 36446,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64811, 0, 3,
                                                                       64481, 36256, 64496,
                                                                       15094, 15112, 36476,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64856, 0, 3,
                                                                       64496, 36266, 64511,
                                                                       15112, 15130, 36506,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64901, 0, 3,
                                                                       64511, 36276, 64526,
                                                                       15130, 15148, 36536,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64946, 0, 3,
                                                                       64526, 36286, 64541,
                                                                       15148, 15166, 36566,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 64991, 0, 3,
                                                                       64541, 36296, 64556,
                                                                       15166, 15184, 36596,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65036, 0, 3,
                                                                       64556, 36306, 64571,
                                                                       15184, 15202, 36626,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65081, 0, 3,
                                                                       64571, 36316, 64586,
                                                                       15202, 15220, 36656,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65126, 0, 3,
                                                                       64586, 36326, 64601,
                                                                       15220, 15238, 36686,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65171, 0, 3,
                                                                       64601, 36336, 64616,
                                                                       15238, 15256, 36716,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65216, 0, 3,
                                                                       64616, 36346, 64631,
                                                                       15256, 15274, 36746,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65261, 0, 3,
                                                                       64631, 36356, 64646,
                                                                       15274, 15292, 36776,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 65306, 0, 3,
                                                                       64646, 36366, 64661,
                                                                       15292, 15310, 36806,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65351, 0, 3,
                                                                       64676, 36386, 64721,
                                                                       15346, 15382, 36836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65441, 0, 3,
                                                                       64721, 36416, 64766,
                                                                       15382, 15418, 36896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65531, 0, 3,
                                                                       64766, 36446, 64811,
                                                                       15418, 15454, 36956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65621, 0, 3,
                                                                       64811, 36476, 64856,
                                                                       15454, 15490, 37016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65711, 0, 3,
                                                                       64856, 36506, 64901,
                                                                       15490, 15526, 37076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65801, 0, 3,
                                                                       64901, 36536, 64946,
                                                                       15526, 15562, 37136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65891, 0, 3,
                                                                       64946, 36566, 64991,
                                                                       15562, 15598, 37196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 65981, 0, 3,
                                                                       64991, 36596, 65036,
                                                                       15598, 15634, 37256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66071, 0, 3,
                                                                       65036, 36626, 65081,
                                                                       15634, 15670, 37316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66161, 0, 3,
                                                                       65081, 36656, 65126,
                                                                       15670, 15706, 37376,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66251, 0, 3,
                                                                       65126, 36686, 65171,
                                                                       15706, 15742, 37436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66341, 0, 3,
                                                                       65171, 36716, 65216,
                                                                       15742, 15778, 37496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66431, 0, 3,
                                                                       65216, 36746, 65261,
                                                                       15778, 15814, 37556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 66521, 0, 3,
                                                                       65261, 36776, 65306,
                                                                       15814, 15850, 37616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 66611, 0, 3,
                                                                       65351, 36836, 65441,
                                                                       15922, 15982, 37676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 66761, 0, 3,
                                                                       65441, 36896, 65531,
                                                                       15982, 16042, 37776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 66911, 0, 3,
                                                                       65531, 36956, 65621,
                                                                       16042, 16102, 37876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67061, 0, 3,
                                                                       65621, 37016, 65711,
                                                                       16102, 16162, 37976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67211, 0, 3,
                                                                       65711, 37076, 65801,
                                                                       16162, 16222, 38076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67361, 0, 3,
                                                                       65801, 37136, 65891,
                                                                       16222, 16282, 38176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67511, 0, 3,
                                                                       65891, 37196, 65981,
                                                                       16282, 16342, 38276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67661, 0, 3,
                                                                       65981, 37256, 66071,
                                                                       16342, 16402, 38376,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67811, 0, 3,
                                                                       66071, 37316, 66161,
                                                                       16402, 16462, 38476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 67961, 0, 3,
                                                                       66161, 37376, 66251,
                                                                       16462, 16522, 38576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 68111, 0, 3,
                                                                       66251, 37436, 66341,
                                                                       16522, 16582, 38676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 68261, 0, 3,
                                                                       66341, 37496, 66431,
                                                                       16582, 16642, 38776,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 68411, 0, 3,
                                                                       66431, 37556, 66521,
                                                                       16642, 16702, 38876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 68561, 0, 3,
                                                                       66611, 37676, 66761,
                                                                       16822, 16912, 38976,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 68786, 0, 3,
                                                                       66761, 37776, 66911,
                                                                       16912, 17002, 39126,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 69011, 0, 3,
                                                                       66911, 37876, 67061,
                                                                       17002, 17092, 39276,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 69236, 0, 3,
                                                                       67061, 37976, 67211,
                                                                       17092, 17182, 39426,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 69461, 0, 3,
                                                                       67211, 38076, 67361,
                                                                       17182, 17272, 39576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 69686, 0, 3,
                                                                       67361, 38176, 67511,
                                                                       17272, 17362, 39726,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 69911, 0, 3,
                                                                       67511, 38276, 67661,
                                                                       17362, 17452, 39876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 70136, 0, 3,
                                                                       67661, 38376, 67811,
                                                                       17452, 17542, 40026,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 70361, 0, 3,
                                                                       67811, 38476, 67961,
                                                                       17542, 17632, 40176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 70586, 0, 3,
                                                                       67961, 38576, 68111,
                                                                       17632, 17722, 40326,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 70811, 0, 3,
                                                                       68111, 38676, 68261,
                                                                       17722, 17812, 40476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 71036, 0, 3,
                                                                       68261, 38776, 68411,
                                                                       17812, 17902, 40626,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 71261, 0, 3,
                                                                       68561, 38976, 68786,
                                                                       18082, 18208, 40776,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 71576, 0, 3,
                                                                       68786, 39126, 69011,
                                                                       18208, 18334, 40986,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 71891, 0, 3,
                                                                       69011, 39276, 69236,
                                                                       18334, 18460, 41196,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 72206, 0, 3,
                                                                       69236, 39426, 69461,
                                                                       18460, 18586, 41406,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 72521, 0, 3,
                                                                       69461, 39576, 69686,
                                                                       18586, 18712, 41616,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 72836, 0, 3,
                                                                       69686, 39726, 69911,
                                                                       18712, 18838, 41826,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 73151, 0, 3,
                                                                       69911, 39876, 70136,
                                                                       18838, 18964, 42036,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 73466, 0, 3,
                                                                       70136, 40026, 70361,
                                                                       18964, 19090, 42246,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 73781, 0, 3,
                                                                       70361, 40176, 70586,
                                                                       19090, 19216, 42456,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 74096, 0, 3,
                                                                       70586, 40326, 70811,
                                                                       19216, 19342, 42666,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 74411, 0, 3,
                                                                       70811, 40476, 71036,
                                                                       19342, 19468, 42876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 74726, 0, 3,
                                                                       71261, 40776, 71576,
                                                                       19720, 19888, 43086,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 75146, 0, 3,
                                                                       71576, 40986, 71891,
                                                                       19888, 20056, 43366,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 75566, 0, 3,
                                                                       71891, 41196, 72206,
                                                                       20056, 20224, 43646,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 75986, 0, 3,
                                                                       72206, 41406, 72521,
                                                                       20224, 20392, 43926,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 76406, 0, 3,
                                                                       72521, 41616, 72836,
                                                                       20392, 20560, 44206,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 76826, 0, 3,
                                                                       72836, 41826, 73151,
                                                                       20560, 20728, 44486,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 77246, 0, 3,
                                                                       73151, 42036, 73466,
                                                                       20728, 20896, 44766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 77666, 0, 3,
                                                                       73466, 42246, 73781,
                                                                       20896, 21064, 45046,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 78086, 0, 3,
                                                                       73781, 42456, 74096,
                                                                       21064, 21232, 45326,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 78506, 0, 3,
                                                                       74096, 42666, 74411,
                                                                       21232, 21400, 45606,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 78926, 0, 3,
                                                                       74726, 43086, 75146,
                                                                       21736, 21952, 45886,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 79466, 0, 3,
                                                                       75146, 43366, 75566,
                                                                       21952, 22168, 46246,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 80006, 0, 3,
                                                                       75566, 43646, 75986,
                                                                       22168, 22384, 46606,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 80546, 0, 3,
                                                                       75986, 43926, 76406,
                                                                       22384, 22600, 46966,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 81086, 0, 3,
                                                                       76406, 44206, 76826,
                                                                       22600, 22816, 47326,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 81626, 0, 3,
                                                                       76826, 44486, 77246,
                                                                       22816, 23032, 47686,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 82166, 0, 3,
                                                                       77246, 44766, 77666,
                                                                       23032, 23248, 48046,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 82706, 0, 3,
                                                                       77666, 45046, 78086,
                                                                       23248, 23464, 48406,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 83246, 0, 3,
                                                                       78086, 45326, 78506,
                                                                       23464, 23680, 48766,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 83786, 0, 3,
                                                                       78926, 45886, 79466,
                                                                       24112, 24382, 49126,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 84461, 0, 3,
                                                                       79466, 46246, 80006,
                                                                       24382, 24652, 49576,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 85136, 0, 3,
                                                                       80006, 46606, 80546,
                                                                       24652, 24922, 50026,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 85811, 0, 3,
                                                                       80546, 46966, 81086,
                                                                       24922, 25192, 50476,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 86486, 0, 3,
                                                                       81086, 47326, 81626,
                                                                       25192, 25462, 50926,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 87161, 0, 3,
                                                                       81626, 47686, 82166,
                                                                       25462, 25732, 51376,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 87836, 0, 3,
                                                                       82166, 48046, 82706,
                                                                       25732, 26002, 51826,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 88511, 0, 3,
                                                                       82706, 48406, 83246,
                                                                       26002, 26272, 52276,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 89186, 0, 3,
                                                                       83786, 49126, 84461,
                                                                       26812, 27142, 52726,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 90011, 0, 3,
                                                                       84461, 49576, 85136,
                                                                       27142, 27472, 53276,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 90836, 0, 3,
                                                                       85136, 50026, 85811,
                                                                       27472, 27802, 53826,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 91661, 0, 3,
                                                                       85811, 50476, 86486,
                                                                       27802, 28132, 54376,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 92486, 0, 3,
                                                                       86486, 50926, 87161,
                                                                       28132, 28462, 54926,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 93311, 0, 3,
                                                                       87161, 51376, 87836,
                                                                       28462, 28792, 55476,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 94136, 0, 3,
                                                                       87836, 51826, 88511,
                                                                       28792, 29122, 56026,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 94961, 0, 3,
                                                                       89186, 52726, 90011,
                                                                       29782, 30178, 56576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 95951, 0, 3,
                                                                       90011, 53276, 90836,
                                                                       30178, 30574, 57236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 96941, 0, 3,
                                                                       90836, 53826, 91661,
                                                                       30574, 30970, 57896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 97931, 0, 3,
                                                                       91661, 54376, 92486,
                                                                       30970, 31366, 58556,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 98921, 0, 3,
                                                                       92486, 54926, 93311,
                                                                       31366, 31762, 59216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 99911, 0, 3,
                                                                       93311, 55476, 94136,
                                                                       31762, 32158, 59876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 100901, 0, 3,
                                                                       94961, 56576, 95951,
                                                                       32950, 33418, 60536,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 102071, 0, 3,
                                                                       95951, 57236, 96941,
                                                                       33418, 33886, 61316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 103241, 0, 3,
                                                                       96941, 57896, 97931,
                                                                       33886, 34354, 62096,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 104411, 0, 3,
                                                                       97931, 58556, 98921,
                                                                       34354, 34822, 62876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 105581, 0, 3,
                                                                       98921, 59216, 99911,
                                                                       34822, 35290, 63656,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106751, 3, 36226,
                                                                       36236, 64466, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106772, 3, 36236,
                                                                       36246, 64481, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106793, 3, 36246,
                                                                       36256, 64496, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106814, 3, 36256,
                                                                       36266, 64511, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106835, 3, 36266,
                                                                       36276, 64526, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106856, 3, 36276,
                                                                       36286, 64541, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106877, 3, 36286,
                                                                       36296, 64556, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106898, 3, 36296,
                                                                       36306, 64571, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106919, 3, 36306,
                                                                       36316, 64586, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106940, 3, 36316,
                                                                       36326, 64601, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106961, 3, 36326,
                                                                       36336, 64616, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 106982, 3, 36336,
                                                                       36346, 64631, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 107003, 3, 36346,
                                                                       36356, 64646, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 107024, 3, 36356,
                                                                       36366, 64661, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107045, 0, 3,
                                                                       106751, 64466, 106772,
                                                                       36386, 36416, 64766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107108, 0, 3,
                                                                       106772, 64481, 106793,
                                                                       36416, 36446, 64811,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107171, 0, 3,
                                                                       106793, 64496, 106814,
                                                                       36446, 36476, 64856,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107234, 0, 3,
                                                                       106814, 64511, 106835,
                                                                       36476, 36506, 64901,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107297, 0, 3,
                                                                       106835, 64526, 106856,
                                                                       36506, 36536, 64946,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107360, 0, 3,
                                                                       106856, 64541, 106877,
                                                                       36536, 36566, 64991,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107423, 0, 3,
                                                                       106877, 64556, 106898,
                                                                       36566, 36596, 65036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107486, 0, 3,
                                                                       106898, 64571, 106919,
                                                                       36596, 36626, 65081,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107549, 0, 3,
                                                                       106919, 64586, 106940,
                                                                       36626, 36656, 65126,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107612, 0, 3,
                                                                       106940, 64601, 106961,
                                                                       36656, 36686, 65171,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107675, 0, 3,
                                                                       106961, 64616, 106982,
                                                                       36686, 36716, 65216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107738, 0, 3,
                                                                       106982, 64631, 107003,
                                                                       36716, 36746, 65261,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 107801, 0, 3,
                                                                       107003, 64646, 107024,
                                                                       36746, 36776, 65306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 107864, 0, 3,
                                                                       107045, 64766, 107108,
                                                                       36836, 36896, 65531,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 107990, 0, 3,
                                                                       107108, 64811, 107171,
                                                                       36896, 36956, 65621,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108116, 0, 3,
                                                                       107171, 64856, 107234,
                                                                       36956, 37016, 65711,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108242, 0, 3,
                                                                       107234, 64901, 107297,
                                                                       37016, 37076, 65801,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108368, 0, 3,
                                                                       107297, 64946, 107360,
                                                                       37076, 37136, 65891,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108494, 0, 3,
                                                                       107360, 64991, 107423,
                                                                       37136, 37196, 65981,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108620, 0, 3,
                                                                       107423, 65036, 107486,
                                                                       37196, 37256, 66071,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108746, 0, 3,
                                                                       107486, 65081, 107549,
                                                                       37256, 37316, 66161,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108872, 0, 3,
                                                                       107549, 65126, 107612,
                                                                       37316, 37376, 66251,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 108998, 0, 3,
                                                                       107612, 65171, 107675,
                                                                       37376, 37436, 66341,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 109124, 0, 3,
                                                                       107675, 65216, 107738,
                                                                       37436, 37496, 66431,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 109250, 0, 3,
                                                                       107738, 65261, 107801,
                                                                       37496, 37556, 66521,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109376, 0, 3,
                                                                       107864, 65531, 107990,
                                                                       37676, 37776, 66911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109586, 0, 3,
                                                                       107990, 65621, 108116,
                                                                       37776, 37876, 67061,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 109796, 0, 3,
                                                                       108116, 65711, 108242,
                                                                       37876, 37976, 67211,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 110006, 0, 3,
                                                                       108242, 65801, 108368,
                                                                       37976, 38076, 67361,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 110216, 0, 3,
                                                                       108368, 65891, 108494,
                                                                       38076, 38176, 67511,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 110426, 0, 3,
                                                                       108494, 65981, 108620,
                                                                       38176, 38276, 67661,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 110636, 0, 3,
                                                                       108620, 66071, 108746,
                                                                       38276, 38376, 67811,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 110846, 0, 3,
                                                                       108746, 66161, 108872,
                                                                       38376, 38476, 67961,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 111056, 0, 3,
                                                                       108872, 66251, 108998,
                                                                       38476, 38576, 68111,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 111266, 0, 3,
                                                                       108998, 66341, 109124,
                                                                       38576, 38676, 68261,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 111476, 0, 3,
                                                                       109124, 66431, 109250,
                                                                       38676, 38776, 68411,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 111686, 0, 3,
                                                                       109376, 66911, 109586,
                                                                       38976, 39126, 69011,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112001, 0, 3,
                                                                       109586, 67061, 109796,
                                                                       39126, 39276, 69236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112316, 0, 3,
                                                                       109796, 67211, 110006,
                                                                       39276, 39426, 69461,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112631, 0, 3,
                                                                       110006, 67361, 110216,
                                                                       39426, 39576, 69686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 112946, 0, 3,
                                                                       110216, 67511, 110426,
                                                                       39576, 39726, 69911,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 113261, 0, 3,
                                                                       110426, 67661, 110636,
                                                                       39726, 39876, 70136,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 113576, 0, 3,
                                                                       110636, 67811, 110846,
                                                                       39876, 40026, 70361,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 113891, 0, 3,
                                                                       110846, 67961, 111056,
                                                                       40026, 40176, 70586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 114206, 0, 3,
                                                                       111056, 68111, 111266,
                                                                       40176, 40326, 70811,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 114521, 0, 3,
                                                                       111266, 68261, 111476,
                                                                       40326, 40476, 71036,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 114836, 0, 3,
                                                                       111686, 69011, 112001,
                                                                       40776, 40986, 71891,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 115277, 0, 3,
                                                                       112001, 69236, 112316,
                                                                       40986, 41196, 72206,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 115718, 0, 3,
                                                                       112316, 69461, 112631,
                                                                       41196, 41406, 72521,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 116159, 0, 3,
                                                                       112631, 69686, 112946,
                                                                       41406, 41616, 72836,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 116600, 0, 3,
                                                                       112946, 69911, 113261,
                                                                       41616, 41826, 73151,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 117041, 0, 3,
                                                                       113261, 70136, 113576,
                                                                       41826, 42036, 73466,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 117482, 0, 3,
                                                                       113576, 70361, 113891,
                                                                       42036, 42246, 73781,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 117923, 0, 3,
                                                                       113891, 70586, 114206,
                                                                       42246, 42456, 74096,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 118364, 0, 3,
                                                                       114206, 70811, 114521,
                                                                       42456, 42666, 74411,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 118805, 0, 3,
                                                                       114836, 71891, 115277,
                                                                       43086, 43366, 75566,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 119393, 0, 3,
                                                                       115277, 72206, 115718,
                                                                       43366, 43646, 75986,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 119981, 0, 3,
                                                                       115718, 72521, 116159,
                                                                       43646, 43926, 76406,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 120569, 0, 3,
                                                                       116159, 72836, 116600,
                                                                       43926, 44206, 76826,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 121157, 0, 3,
                                                                       116600, 73151, 117041,
                                                                       44206, 44486, 77246,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 121745, 0, 3,
                                                                       117041, 73466, 117482,
                                                                       44486, 44766, 77666,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 122333, 0, 3,
                                                                       117482, 73781, 117923,
                                                                       44766, 45046, 78086,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 122921, 0, 3,
                                                                       117923, 74096, 118364,
                                                                       45046, 45326, 78506,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 123509, 0, 3,
                                                                       118805, 75566, 119393,
                                                                       45886, 46246, 80006,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 124265, 0, 3,
                                                                       119393, 75986, 119981,
                                                                       46246, 46606, 80546,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 125021, 0, 3,
                                                                       119981, 76406, 120569,
                                                                       46606, 46966, 81086,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 125777, 0, 3,
                                                                       120569, 76826, 121157,
                                                                       46966, 47326, 81626,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 126533, 0, 3,
                                                                       121157, 77246, 121745,
                                                                       47326, 47686, 82166,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 127289, 0, 3,
                                                                       121745, 77666, 122333,
                                                                       47686, 48046, 82706,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 128045, 0, 3,
                                                                       122333, 78086, 122921,
                                                                       48046, 48406, 83246,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 128801, 0, 3,
                                                                       123509, 80006, 124265,
                                                                       49126, 49576, 85136,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 129746, 0, 3,
                                                                       124265, 80546, 125021,
                                                                       49576, 50026, 85811,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 130691, 0, 3,
                                                                       125021, 81086, 125777,
                                                                       50026, 50476, 86486,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 131636, 0, 3,
                                                                       125777, 81626, 126533,
                                                                       50476, 50926, 87161,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 132581, 0, 3,
                                                                       126533, 82166, 127289,
                                                                       50926, 51376, 87836,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 133526, 0, 3,
                                                                       127289, 82706, 128045,
                                                                       51376, 51826, 88511,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 134471, 0, 3,
                                                                       128801, 85136, 129746,
                                                                       52726, 53276, 90836,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 135626, 0, 3,
                                                                       129746, 85811, 130691,
                                                                       53276, 53826, 91661,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 136781, 0, 3,
                                                                       130691, 86486, 131636,
                                                                       53826, 54376, 92486,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 137936, 0, 3,
                                                                       131636, 87161, 132581,
                                                                       54376, 54926, 93311,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 139091, 0, 3,
                                                                       132581, 87836, 133526,
                                                                       54926, 55476, 94136,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 140246, 0, 3,
                                                                       134471, 90836, 135626,
                                                                       56576, 57236, 96941,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 141632, 0, 3,
                                                                       135626, 91661, 136781,
                                                                       57236, 57896, 97931,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 143018, 0, 3,
                                                                       136781, 92486, 137936,
                                                                       57896, 58556, 98921,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 144404, 0, 3,
                                                                       137936, 93311, 139091,
                                                                       58556, 59216, 99911,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 145790, 0, 3,
                                                                       140246, 96941, 141632,
                                                                       60536, 61316, 103241,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 147428, 0, 3,
                                                                       141632, 97931, 143018,
                                                                       61316, 62096, 104411,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 149066, 0, 3,
                                                                       143018, 98921, 144404,
                                                                       62096, 62876, 105581,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150704, 3, 64436,
                                                                       64451, 106751, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150732, 3, 64451,
                                                                       64466, 106772, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150760, 3, 64466,
                                                                       64481, 106793, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150788, 3, 64481,
                                                                       64496, 106814, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150816, 3, 64496,
                                                                       64511, 106835, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150844, 3, 64511,
                                                                       64526, 106856, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150872, 3, 64526,
                                                                       64541, 106877, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150900, 3, 64541,
                                                                       64556, 106898, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150928, 3, 64556,
                                                                       64571, 106919, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150956, 3, 64571,
                                                                       64586, 106940, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 150984, 3, 64586,
                                                                       64601, 106961, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 151012, 3, 64601,
                                                                       64616, 106982, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 151040, 3, 64616,
                                                                       64631, 107003, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 151068, 3, 64631,
                                                                       64646, 107024, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151096, 0, 3,
                                                                       150704, 106751, 150732,
                                                                       64676, 64721, 107045,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151180, 0, 3,
                                                                       150732, 106772, 150760,
                                                                       64721, 64766, 107108,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151264, 0, 3,
                                                                       150760, 106793, 150788,
                                                                       64766, 64811, 107171,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151348, 0, 3,
                                                                       150788, 106814, 150816,
                                                                       64811, 64856, 107234,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151432, 0, 3,
                                                                       150816, 106835, 150844,
                                                                       64856, 64901, 107297,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151516, 0, 3,
                                                                       150844, 106856, 150872,
                                                                       64901, 64946, 107360,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151600, 0, 3,
                                                                       150872, 106877, 150900,
                                                                       64946, 64991, 107423,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151684, 0, 3,
                                                                       150900, 106898, 150928,
                                                                       64991, 65036, 107486,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151768, 0, 3,
                                                                       150928, 106919, 150956,
                                                                       65036, 65081, 107549,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151852, 0, 3,
                                                                       150956, 106940, 150984,
                                                                       65081, 65126, 107612,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 151936, 0, 3,
                                                                       150984, 106961, 151012,
                                                                       65126, 65171, 107675,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 152020, 0, 3,
                                                                       151012, 106982, 151040,
                                                                       65171, 65216, 107738,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 152104, 0, 3,
                                                                       151040, 107003, 151068,
                                                                       65216, 65261, 107801,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152188, 0, 3,
                                                                       151096, 107045, 151180,
                                                                       65351, 65441, 107864,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152356, 0, 3,
                                                                       151180, 107108, 151264,
                                                                       65441, 65531, 107990,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152524, 0, 3,
                                                                       151264, 107171, 151348,
                                                                       65531, 65621, 108116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152692, 0, 3,
                                                                       151348, 107234, 151432,
                                                                       65621, 65711, 108242,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 152860, 0, 3,
                                                                       151432, 107297, 151516,
                                                                       65711, 65801, 108368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153028, 0, 3,
                                                                       151516, 107360, 151600,
                                                                       65801, 65891, 108494,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153196, 0, 3,
                                                                       151600, 107423, 151684,
                                                                       65891, 65981, 108620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153364, 0, 3,
                                                                       151684, 107486, 151768,
                                                                       65981, 66071, 108746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153532, 0, 3,
                                                                       151768, 107549, 151852,
                                                                       66071, 66161, 108872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153700, 0, 3,
                                                                       151852, 107612, 151936,
                                                                       66161, 66251, 108998,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 153868, 0, 3,
                                                                       151936, 107675, 152020,
                                                                       66251, 66341, 109124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 154036, 0, 3,
                                                                       152020, 107738, 152104,
                                                                       66341, 66431, 109250,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 154204, 0, 3,
                                                                       152188, 107864, 152356,
                                                                       66611, 66761, 109376,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 154484, 0, 3,
                                                                       152356, 107990, 152524,
                                                                       66761, 66911, 109586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 154764, 0, 3,
                                                                       152524, 108116, 152692,
                                                                       66911, 67061, 109796,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155044, 0, 3,
                                                                       152692, 108242, 152860,
                                                                       67061, 67211, 110006,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155324, 0, 3,
                                                                       152860, 108368, 153028,
                                                                       67211, 67361, 110216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155604, 0, 3,
                                                                       153028, 108494, 153196,
                                                                       67361, 67511, 110426,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 155884, 0, 3,
                                                                       153196, 108620, 153364,
                                                                       67511, 67661, 110636,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 156164, 0, 3,
                                                                       153364, 108746, 153532,
                                                                       67661, 67811, 110846,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 156444, 0, 3,
                                                                       153532, 108872, 153700,
                                                                       67811, 67961, 111056,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 156724, 0, 3,
                                                                       153700, 108998, 153868,
                                                                       67961, 68111, 111266,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 157004, 0, 3,
                                                                       153868, 109124, 154036,
                                                                       68111, 68261, 111476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 157284, 0, 3,
                                                                       154204, 109376, 154484,
                                                                       68561, 68786, 111686,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 157704, 0, 3,
                                                                       154484, 109586, 154764,
                                                                       68786, 69011, 112001,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 158124, 0, 3,
                                                                       154764, 109796, 155044,
                                                                       69011, 69236, 112316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 158544, 0, 3,
                                                                       155044, 110006, 155324,
                                                                       69236, 69461, 112631,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 158964, 0, 3,
                                                                       155324, 110216, 155604,
                                                                       69461, 69686, 112946,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 159384, 0, 3,
                                                                       155604, 110426, 155884,
                                                                       69686, 69911, 113261,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 159804, 0, 3,
                                                                       155884, 110636, 156164,
                                                                       69911, 70136, 113576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 160224, 0, 3,
                                                                       156164, 110846, 156444,
                                                                       70136, 70361, 113891,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 160644, 0, 3,
                                                                       156444, 111056, 156724,
                                                                       70361, 70586, 114206,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 161064, 0, 3,
                                                                       156724, 111266, 157004,
                                                                       70586, 70811, 114521,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 161484, 0, 3,
                                                                       157284, 111686, 157704,
                                                                       71261, 71576, 114836,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 162072, 0, 3,
                                                                       157704, 112001, 158124,
                                                                       71576, 71891, 115277,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 162660, 0, 3,
                                                                       158124, 112316, 158544,
                                                                       71891, 72206, 115718,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 163248, 0, 3,
                                                                       158544, 112631, 158964,
                                                                       72206, 72521, 116159,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 163836, 0, 3,
                                                                       158964, 112946, 159384,
                                                                       72521, 72836, 116600,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 164424, 0, 3,
                                                                       159384, 113261, 159804,
                                                                       72836, 73151, 117041,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 165012, 0, 3,
                                                                       159804, 113576, 160224,
                                                                       73151, 73466, 117482,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 165600, 0, 3,
                                                                       160224, 113891, 160644,
                                                                       73466, 73781, 117923,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 166188, 0, 3,
                                                                       160644, 114206, 161064,
                                                                       73781, 74096, 118364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 166776, 0, 3,
                                                                       161484, 114836, 162072,
                                                                       74726, 75146, 118805,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 167560, 0, 3,
                                                                       162072, 115277, 162660,
                                                                       75146, 75566, 119393,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 168344, 0, 3,
                                                                       162660, 115718, 163248,
                                                                       75566, 75986, 119981,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 169128, 0, 3,
                                                                       163248, 116159, 163836,
                                                                       75986, 76406, 120569,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 169912, 0, 3,
                                                                       163836, 116600, 164424,
                                                                       76406, 76826, 121157,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 170696, 0, 3,
                                                                       164424, 117041, 165012,
                                                                       76826, 77246, 121745,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 171480, 0, 3,
                                                                       165012, 117482, 165600,
                                                                       77246, 77666, 122333,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 172264, 0, 3,
                                                                       165600, 117923, 166188,
                                                                       77666, 78086, 122921,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 173048, 0, 3,
                                                                       166776, 118805, 167560,
                                                                       78926, 79466, 123509,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 174056, 0, 3,
                                                                       167560, 119393, 168344,
                                                                       79466, 80006, 124265,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 175064, 0, 3,
                                                                       168344, 119981, 169128,
                                                                       80006, 80546, 125021,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 176072, 0, 3,
                                                                       169128, 120569, 169912,
                                                                       80546, 81086, 125777,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 177080, 0, 3,
                                                                       169912, 121157, 170696,
                                                                       81086, 81626, 126533,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 178088, 0, 3,
                                                                       170696, 121745, 171480,
                                                                       81626, 82166, 127289,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 179096, 0, 3,
                                                                       171480, 122333, 172264,
                                                                       82166, 82706, 128045,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 180104, 0, 3,
                                                                       173048, 123509, 174056,
                                                                       83786, 84461, 128801,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 181364, 0, 3,
                                                                       174056, 124265, 175064,
                                                                       84461, 85136, 129746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 182624, 0, 3,
                                                                       175064, 125021, 176072,
                                                                       85136, 85811, 130691,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 183884, 0, 3,
                                                                       176072, 125777, 177080,
                                                                       85811, 86486, 131636,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 185144, 0, 3,
                                                                       177080, 126533, 178088,
                                                                       86486, 87161, 132581,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 186404, 0, 3,
                                                                       178088, 127289, 179096,
                                                                       87161, 87836, 133526,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 187664, 0, 3,
                                                                       180104, 128801, 181364,
                                                                       89186, 90011, 134471,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 189204, 0, 3,
                                                                       181364, 129746, 182624,
                                                                       90011, 90836, 135626,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 190744, 0, 3,
                                                                       182624, 130691, 183884,
                                                                       90836, 91661, 136781,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 192284, 0, 3,
                                                                       183884, 131636, 185144,
                                                                       91661, 92486, 137936,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 193824, 0, 3,
                                                                       185144, 132581, 186404,
                                                                       92486, 93311, 139091,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 195364, 0, 3,
                                                                       187664, 134471, 189204,
                                                                       94961, 95951, 140246,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 197212, 0, 3,
                                                                       189204, 135626, 190744,
                                                                       95951, 96941, 141632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 199060, 0, 3,
                                                                       190744, 136781, 192284,
                                                                       96941, 97931, 143018,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 200908, 0, 3,
                                                                       192284, 137936, 193824,
                                                                       97931, 98921, 144404,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 202756, 0, 3,
                                                                       195364, 140246, 197212,
                                                                       100901, 102071, 145790,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 204940, 0, 3,
                                                                       197212, 141632, 199060,
                                                                       102071, 103241, 147428,
                                                                       ncols, gamma, p, q);

                    compute_prim_soi_three_center_electron_repulsion_0(buffer, 207124, 0, 3,
                                                                       199060, 143018, 200908,
                                                                       103241, 104411, 149066,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209308, 3, 106751,
                                                                       106772, 150760, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209344, 3, 106772,
                                                                       106793, 150788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209380, 3, 106793,
                                                                       106814, 150816, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209416, 3, 106814,
                                                                       106835, 150844, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209452, 3, 106835,
                                                                       106856, 150872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209488, 3, 106856,
                                                                       106877, 150900, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209524, 3, 106877,
                                                                       106898, 150928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209560, 3, 106898,
                                                                       106919, 150956, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209596, 3, 106919,
                                                                       106940, 150984, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209632, 3, 106940,
                                                                       106961, 151012, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209668, 3, 106961,
                                                                       106982, 151040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 209704, 3, 106982,
                                                                       107003, 151068, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 209740, 0, 3,
                                                                       209308, 150760, 209344,
                                                                       107045, 107108, 151264,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 209848, 0, 3,
                                                                       209344, 150788, 209380,
                                                                       107108, 107171, 151348,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 209956, 0, 3,
                                                                       209380, 150816, 209416,
                                                                       107171, 107234, 151432,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210064, 0, 3,
                                                                       209416, 150844, 209452,
                                                                       107234, 107297, 151516,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210172, 0, 3,
                                                                       209452, 150872, 209488,
                                                                       107297, 107360, 151600,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210280, 0, 3,
                                                                       209488, 150900, 209524,
                                                                       107360, 107423, 151684,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210388, 0, 3,
                                                                       209524, 150928, 209560,
                                                                       107423, 107486, 151768,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210496, 0, 3,
                                                                       209560, 150956, 209596,
                                                                       107486, 107549, 151852,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210604, 0, 3,
                                                                       209596, 150984, 209632,
                                                                       107549, 107612, 151936,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210712, 0, 3,
                                                                       209632, 151012, 209668,
                                                                       107612, 107675, 152020,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 210820, 0, 3,
                                                                       209668, 151040, 209704,
                                                                       107675, 107738, 152104,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 210928, 0, 3,
                                                                       209740, 151264, 209848,
                                                                       107864, 107990, 152524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 211144, 0, 3,
                                                                       209848, 151348, 209956,
                                                                       107990, 108116, 152692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 211360, 0, 3,
                                                                       209956, 151432, 210064,
                                                                       108116, 108242, 152860,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 211576, 0, 3,
                                                                       210064, 151516, 210172,
                                                                       108242, 108368, 153028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 211792, 0, 3,
                                                                       210172, 151600, 210280,
                                                                       108368, 108494, 153196,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 212008, 0, 3,
                                                                       210280, 151684, 210388,
                                                                       108494, 108620, 153364,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 212224, 0, 3,
                                                                       210388, 151768, 210496,
                                                                       108620, 108746, 153532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 212440, 0, 3,
                                                                       210496, 151852, 210604,
                                                                       108746, 108872, 153700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 212656, 0, 3,
                                                                       210604, 151936, 210712,
                                                                       108872, 108998, 153868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 212872, 0, 3,
                                                                       210712, 152020, 210820,
                                                                       108998, 109124, 154036,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 213088, 0, 3,
                                                                       210928, 152524, 211144,
                                                                       109376, 109586, 154764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 213448, 0, 3,
                                                                       211144, 152692, 211360,
                                                                       109586, 109796, 155044,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 213808, 0, 3,
                                                                       211360, 152860, 211576,
                                                                       109796, 110006, 155324,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 214168, 0, 3,
                                                                       211576, 153028, 211792,
                                                                       110006, 110216, 155604,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 214528, 0, 3,
                                                                       211792, 153196, 212008,
                                                                       110216, 110426, 155884,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 214888, 0, 3,
                                                                       212008, 153364, 212224,
                                                                       110426, 110636, 156164,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 215248, 0, 3,
                                                                       212224, 153532, 212440,
                                                                       110636, 110846, 156444,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 215608, 0, 3,
                                                                       212440, 153700, 212656,
                                                                       110846, 111056, 156724,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 215968, 0, 3,
                                                                       212656, 153868, 212872,
                                                                       111056, 111266, 157004,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 216328, 0, 3,
                                                                       213088, 154764, 213448,
                                                                       111686, 112001, 158124,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 216868, 0, 3,
                                                                       213448, 155044, 213808,
                                                                       112001, 112316, 158544,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 217408, 0, 3,
                                                                       213808, 155324, 214168,
                                                                       112316, 112631, 158964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 217948, 0, 3,
                                                                       214168, 155604, 214528,
                                                                       112631, 112946, 159384,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 218488, 0, 3,
                                                                       214528, 155884, 214888,
                                                                       112946, 113261, 159804,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 219028, 0, 3,
                                                                       214888, 156164, 215248,
                                                                       113261, 113576, 160224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 219568, 0, 3,
                                                                       215248, 156444, 215608,
                                                                       113576, 113891, 160644,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 220108, 0, 3,
                                                                       215608, 156724, 215968,
                                                                       113891, 114206, 161064,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 220648, 0, 3,
                                                                       216328, 158124, 216868,
                                                                       114836, 115277, 162660,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 221404, 0, 3,
                                                                       216868, 158544, 217408,
                                                                       115277, 115718, 163248,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 222160, 0, 3,
                                                                       217408, 158964, 217948,
                                                                       115718, 116159, 163836,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 222916, 0, 3,
                                                                       217948, 159384, 218488,
                                                                       116159, 116600, 164424,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 223672, 0, 3,
                                                                       218488, 159804, 219028,
                                                                       116600, 117041, 165012,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 224428, 0, 3,
                                                                       219028, 160224, 219568,
                                                                       117041, 117482, 165600,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 225184, 0, 3,
                                                                       219568, 160644, 220108,
                                                                       117482, 117923, 166188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 225940, 0, 3,
                                                                       220648, 162660, 221404,
                                                                       118805, 119393, 168344,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 226948, 0, 3,
                                                                       221404, 163248, 222160,
                                                                       119393, 119981, 169128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 227956, 0, 3,
                                                                       222160, 163836, 222916,
                                                                       119981, 120569, 169912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 228964, 0, 3,
                                                                       222916, 164424, 223672,
                                                                       120569, 121157, 170696,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 229972, 0, 3,
                                                                       223672, 165012, 224428,
                                                                       121157, 121745, 171480,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 230980, 0, 3,
                                                                       224428, 165600, 225184,
                                                                       121745, 122333, 172264,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 231988, 0, 3,
                                                                       225940, 168344, 226948,
                                                                       123509, 124265, 175064,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 233284, 0, 3,
                                                                       226948, 169128, 227956,
                                                                       124265, 125021, 176072,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 234580, 0, 3,
                                                                       227956, 169912, 228964,
                                                                       125021, 125777, 177080,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 235876, 0, 3,
                                                                       228964, 170696, 229972,
                                                                       125777, 126533, 178088,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 237172, 0, 3,
                                                                       229972, 171480, 230980,
                                                                       126533, 127289, 179096,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 238468, 0, 3,
                                                                       231988, 175064, 233284,
                                                                       128801, 129746, 182624,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 240088, 0, 3,
                                                                       233284, 176072, 234580,
                                                                       129746, 130691, 183884,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 241708, 0, 3,
                                                                       234580, 177080, 235876,
                                                                       130691, 131636, 185144,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 243328, 0, 3,
                                                                       235876, 178088, 237172,
                                                                       131636, 132581, 186404,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 244948, 0, 3,
                                                                       238468, 182624, 240088,
                                                                       134471, 135626, 190744,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 246928, 0, 3,
                                                                       240088, 183884, 241708,
                                                                       135626, 136781, 192284,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 248908, 0, 3,
                                                                       241708, 185144, 243328,
                                                                       136781, 137936, 193824,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 250888, 0, 3,
                                                                       244948, 190744, 246928,
                                                                       140246, 141632, 199060,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 253264, 0, 3,
                                                                       246928, 192284, 248908,
                                                                       141632, 143018, 200908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sok_three_center_electron_repulsion_0(buffer, 255640, 0, 3,
                                                                       250888, 199060, 253264,
                                                                       145790, 147428, 207124,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258448, 3, 150704,
                                                                       150732, 209308, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258493, 3, 150732,
                                                                       150760, 209344, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258538, 3, 150760,
                                                                       150788, 209380, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258583, 3, 150788,
                                                                       150816, 209416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258628, 3, 150816,
                                                                       150844, 209452, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258673, 3, 150844,
                                                                       150872, 209488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258718, 3, 150872,
                                                                       150900, 209524, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258763, 3, 150900,
                                                                       150928, 209560, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258808, 3, 150928,
                                                                       150956, 209596, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258853, 3, 150956,
                                                                       150984, 209632, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258898, 3, 150984,
                                                                       151012, 209668, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 258943, 3, 151012,
                                                                       151040, 209704, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 258988, 0, 3,
                                                                       258448, 209308, 258493,
                                                                       151096, 151180, 209740,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259123, 0, 3,
                                                                       258493, 209344, 258538,
                                                                       151180, 151264, 209848,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259258, 0, 3,
                                                                       258538, 209380, 258583,
                                                                       151264, 151348, 209956,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259393, 0, 3,
                                                                       258583, 209416, 258628,
                                                                       151348, 151432, 210064,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259528, 0, 3,
                                                                       258628, 209452, 258673,
                                                                       151432, 151516, 210172,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259663, 0, 3,
                                                                       258673, 209488, 258718,
                                                                       151516, 151600, 210280,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259798, 0, 3,
                                                                       258718, 209524, 258763,
                                                                       151600, 151684, 210388,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 259933, 0, 3,
                                                                       258763, 209560, 258808,
                                                                       151684, 151768, 210496,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 260068, 0, 3,
                                                                       258808, 209596, 258853,
                                                                       151768, 151852, 210604,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 260203, 0, 3,
                                                                       258853, 209632, 258898,
                                                                       151852, 151936, 210712,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 260338, 0, 3,
                                                                       258898, 209668, 258943,
                                                                       151936, 152020, 210820,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 260473, 0, 3,
                                                                       258988, 209740, 259123,
                                                                       152188, 152356, 210928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 260743, 0, 3,
                                                                       259123, 209848, 259258,
                                                                       152356, 152524, 211144,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 261013, 0, 3,
                                                                       259258, 209956, 259393,
                                                                       152524, 152692, 211360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 261283, 0, 3,
                                                                       259393, 210064, 259528,
                                                                       152692, 152860, 211576,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 261553, 0, 3,
                                                                       259528, 210172, 259663,
                                                                       152860, 153028, 211792,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 261823, 0, 3,
                                                                       259663, 210280, 259798,
                                                                       153028, 153196, 212008,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 262093, 0, 3,
                                                                       259798, 210388, 259933,
                                                                       153196, 153364, 212224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 262363, 0, 3,
                                                                       259933, 210496, 260068,
                                                                       153364, 153532, 212440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 262633, 0, 3,
                                                                       260068, 210604, 260203,
                                                                       153532, 153700, 212656,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 262903, 0, 3,
                                                                       260203, 210712, 260338,
                                                                       153700, 153868, 212872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 263173, 0, 3,
                                                                       260473, 210928, 260743,
                                                                       154204, 154484, 213088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 263623, 0, 3,
                                                                       260743, 211144, 261013,
                                                                       154484, 154764, 213448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 264073, 0, 3,
                                                                       261013, 211360, 261283,
                                                                       154764, 155044, 213808,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 264523, 0, 3,
                                                                       261283, 211576, 261553,
                                                                       155044, 155324, 214168,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 264973, 0, 3,
                                                                       261553, 211792, 261823,
                                                                       155324, 155604, 214528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 265423, 0, 3,
                                                                       261823, 212008, 262093,
                                                                       155604, 155884, 214888,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 265873, 0, 3,
                                                                       262093, 212224, 262363,
                                                                       155884, 156164, 215248,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 266323, 0, 3,
                                                                       262363, 212440, 262633,
                                                                       156164, 156444, 215608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 266773, 0, 3,
                                                                       262633, 212656, 262903,
                                                                       156444, 156724, 215968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 267223, 0, 3,
                                                                       263173, 213088, 263623,
                                                                       157284, 157704, 216328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 267898, 0, 3,
                                                                       263623, 213448, 264073,
                                                                       157704, 158124, 216868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 268573, 0, 3,
                                                                       264073, 213808, 264523,
                                                                       158124, 158544, 217408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 269248, 0, 3,
                                                                       264523, 214168, 264973,
                                                                       158544, 158964, 217948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 269923, 0, 3,
                                                                       264973, 214528, 265423,
                                                                       158964, 159384, 218488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 270598, 0, 3,
                                                                       265423, 214888, 265873,
                                                                       159384, 159804, 219028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 271273, 0, 3,
                                                                       265873, 215248, 266323,
                                                                       159804, 160224, 219568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 271948, 0, 3,
                                                                       266323, 215608, 266773,
                                                                       160224, 160644, 220108,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 272623, 0, 3,
                                                                       267223, 216328, 267898,
                                                                       161484, 162072, 220648,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 273568, 0, 3,
                                                                       267898, 216868, 268573,
                                                                       162072, 162660, 221404,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 274513, 0, 3,
                                                                       268573, 217408, 269248,
                                                                       162660, 163248, 222160,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 275458, 0, 3,
                                                                       269248, 217948, 269923,
                                                                       163248, 163836, 222916,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 276403, 0, 3,
                                                                       269923, 218488, 270598,
                                                                       163836, 164424, 223672,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 277348, 0, 3,
                                                                       270598, 219028, 271273,
                                                                       164424, 165012, 224428,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 278293, 0, 3,
                                                                       271273, 219568, 271948,
                                                                       165012, 165600, 225184,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 279238, 0, 3,
                                                                       272623, 220648, 273568,
                                                                       166776, 167560, 225940,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 280498, 0, 3,
                                                                       273568, 221404, 274513,
                                                                       167560, 168344, 226948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 281758, 0, 3,
                                                                       274513, 222160, 275458,
                                                                       168344, 169128, 227956,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 283018, 0, 3,
                                                                       275458, 222916, 276403,
                                                                       169128, 169912, 228964,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 284278, 0, 3,
                                                                       276403, 223672, 277348,
                                                                       169912, 170696, 229972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 285538, 0, 3,
                                                                       277348, 224428, 278293,
                                                                       170696, 171480, 230980,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 286798, 0, 3,
                                                                       279238, 225940, 280498,
                                                                       173048, 174056, 231988,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 288418, 0, 3,
                                                                       280498, 226948, 281758,
                                                                       174056, 175064, 233284,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 290038, 0, 3,
                                                                       281758, 227956, 283018,
                                                                       175064, 176072, 234580,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 291658, 0, 3,
                                                                       283018, 228964, 284278,
                                                                       176072, 177080, 235876,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 293278, 0, 3,
                                                                       284278, 229972, 285538,
                                                                       177080, 178088, 237172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 294898, 0, 3,
                                                                       286798, 231988, 288418,
                                                                       180104, 181364, 238468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 296923, 0, 3,
                                                                       288418, 233284, 290038,
                                                                       181364, 182624, 240088,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 298948, 0, 3,
                                                                       290038, 234580, 291658,
                                                                       182624, 183884, 241708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sll_three_center_electron_repulsion_0(buffer, 300973, 0, 3,
                                                                       291658, 235876, 293278,
                                                                       183884, 185144, 243328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 302998, 0, 3,
                                                                       294898, 238468, 296923,
                                                                       187664, 189204, 244948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 305473, 0, 3,
                                                                       296923, 240088, 298948,
                                                                       189204, 190744, 246928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sml_three_center_electron_repulsion_0(buffer, 307948, 0, 3,
                                                                       298948, 241708, 300973,
                                                                       190744, 192284, 248908,
                                                                       ncols, gamma, p, q);

                    compute_prim_snl_three_center_electron_repulsion_0(buffer, 310423, 0, 3,
                                                                       302998, 244948, 305473,
                                                                       195364, 197212, 250888,
                                                                       ncols, gamma, p, q);

                    compute_prim_snl_three_center_electron_repulsion_0(buffer, 313393, 0, 3,
                                                                       305473, 246928, 307948,
                                                                       197212, 199060, 253264,
                                                                       ncols, gamma, p, q);

                    compute_prim_sol_three_center_electron_repulsion_0(buffer, 316363, 0, 3,
                                                                       310423, 250888, 313393,
                                                                       202756, 204940, 255640,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 319873, 279238, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 321609, 286798, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 323841, 294898, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 326631, 302998, 2475, ncols);

                    simdfunc::contract_primitives(buffer, 330041, 310423, 2970, ncols);

                    simdfunc::contract_primitives(buffer, 334133, 316363, 3510, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 321133, 319873, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 323229, 321609, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 325866, 323841, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 329106, 326631, 55, 1, nmax);

        simdtrf::transform_l_inner(buffer, 333011, 330041, 66, 1, nmax);

        simdtrf::transform_l_inner(buffer, 337643, 334133, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 338969, 321133, 323229, 17, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 340397, 323229, 325866, 17, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 342233, 325866, 329106, 17, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 344528, 329106, 333011, 17, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 347333, 333011, 337643, 17, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 350699, 338969, 340397, 17, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 353555, 340397, 342233, 17, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 357227, 342233, 344528, 17, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 361817, 344528, 347333, 17, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 367427, 350699, 353555, 17, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 372187, 353555, 357227, 17, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 378307, 357227, 361817, 17, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 385957, 367427, 372187, 17, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 393097, 372187, 378307, 17, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 402277, 385957, 393097, 17, nmax);

        simdtrf::transform_i_inner(buffer, 412273, 402277, 21, 17, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 412273, 221, nmax);
    }

    for (size_t m = 0; m < 2431; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
