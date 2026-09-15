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


#include "SimdThreeCenterElectronRepulsionRecGIK.hpp"

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
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_gik_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_gik_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 213615, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1755 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 213615, 165870, 10740, dimensions);

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
                                                        16, 17}, ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 73, 0, 3, 8, 9,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 79, 0, 3, 9, 10,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 85, 0, 3, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 91, 0, 3, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 97, 0, 3, 12, 13,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 103, 0, 3, 13, 14,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 109, 0, 3, 14, 15,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 115, 0, 3, 15, 16,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 121, 0, 3, 16, 17,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 127, 0, 3, 17, 18,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 133, 0, 3, 18, 19,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 139, 0, 3, 19, 20,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 145, 0, 3, 20, 21,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 151, 0, 3, 21, 22,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 157, 0, 3, 22, 23,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 25, 28,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 28, 31,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 31, 34,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 34, 37,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 37, 40,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 40, 43,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 43, 46,
                                                                       109, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 233, 0, 3, 46, 49,
                                                                       115, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 49, 52,
                                                                       121, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 52, 55,
                                                                       127, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 263, 0, 3, 55, 58,
                                                                       133, 139, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 58, 61,
                                                                       139, 145, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 283, 0, 3, 61, 64,
                                                                       145, 151, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 293, 0, 3, 64, 67,
                                                                       151, 157, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 303, 0, 3, 73, 79,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 79, 85,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 333, 0, 3, 85, 91,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 91, 97,
                                                                       193, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 363, 0, 3, 97,
                                                                       103, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 103,
                                                                       109, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 393, 0, 3, 109,
                                                                       115, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 408, 0, 3, 115,
                                                                       121, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 423, 0, 3, 121,
                                                                       127, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 438, 0, 3, 127,
                                                                       133, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 453, 0, 3, 133,
                                                                       139, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 468, 0, 3, 139,
                                                                       145, 273, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 483, 0, 3, 145,
                                                                       151, 283, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 163,
                                                                       173, 303, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 519, 0, 3, 173,
                                                                       183, 318, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 540, 0, 3, 183,
                                                                       193, 333, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 561, 0, 3, 193,
                                                                       203, 348, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 582, 0, 3, 203,
                                                                       213, 363, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 603, 0, 3, 213,
                                                                       223, 378, 393, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 624, 0, 3, 223,
                                                                       233, 393, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 645, 0, 3, 233,
                                                                       243, 408, 423, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 666, 0, 3, 243,
                                                                       253, 423, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 687, 0, 3, 253,
                                                                       263, 438, 453, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 708, 0, 3, 263,
                                                                       273, 453, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 729, 0, 3, 273,
                                                                       283, 468, 483, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 750, 0, 3, 303,
                                                                       318, 498, 519, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 778, 0, 3, 318,
                                                                       333, 519, 540, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 806, 0, 3, 333,
                                                                       348, 540, 561, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 834, 0, 3, 348,
                                                                       363, 561, 582, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 862, 0, 3, 363,
                                                                       378, 582, 603, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 890, 0, 3, 378,
                                                                       393, 603, 624, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 918, 0, 3, 393,
                                                                       408, 624, 645, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 946, 0, 3, 408,
                                                                       423, 645, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 974, 0, 3, 423,
                                                                       438, 666, 687, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 438,
                                                                       453, 687, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 453,
                                                                       468, 708, 729, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 498,
                                                                       519, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1094, 0, 3, 519,
                                                                       540, 778, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1130, 0, 3, 540,
                                                                       561, 806, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1166, 0, 3, 561,
                                                                       582, 834, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1202, 0, 3, 582,
                                                                       603, 862, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1238, 0, 3, 603,
                                                                       624, 890, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1274, 0, 3, 624,
                                                                       645, 918, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 645,
                                                                       666, 946, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1346, 0, 3, 666,
                                                                       687, 974, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1382, 0, 3, 687,
                                                                       708, 1002, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1418, 0, 3, 750,
                                                                       778, 1058, 1094, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1463, 0, 3, 778,
                                                                       806, 1094, 1130, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1508, 0, 3, 806,
                                                                       834, 1130, 1166, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1553, 0, 3, 834,
                                                                       862, 1166, 1202, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1598, 0, 3, 862,
                                                                       890, 1202, 1238, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1643, 0, 3, 890,
                                                                       918, 1238, 1274, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1688, 0, 3, 918,
                                                                       946, 1274, 1310, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1733, 0, 3, 946,
                                                                       974, 1310, 1346, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1778, 0, 3, 974,
                                                                       1002, 1346, 1382, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1823, 0, 3, 1058,
                                                                       1094, 1418, 1463, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1878, 0, 3, 1094,
                                                                       1130, 1463, 1508, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1933, 0, 3, 1130,
                                                                       1166, 1508, 1553, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1166,
                                                                       1202, 1553, 1598, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2043, 0, 3, 1202,
                                                                       1238, 1598, 1643, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2098, 0, 3, 1238,
                                                                       1274, 1643, 1688, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1274,
                                                                       1310, 1688, 1733, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 2208, 0, 3, 1310,
                                                                       1346, 1733, 1778, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2263, 0, 3, 1418,
                                                                       1463, 1823, 1878, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2329, 0, 3, 1463,
                                                                       1508, 1878, 1933, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2395, 0, 3, 1508,
                                                                       1553, 1933, 1988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2461, 0, 3, 1553,
                                                                       1598, 1988, 2043, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2527, 0, 3, 1598,
                                                                       1643, 2043, 2098, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2593, 0, 3, 1643,
                                                                       1688, 2098, 2153, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2659, 0, 3, 1688,
                                                                       1733, 2153, 2208, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2725, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2728, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2731, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2734, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2737, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2740, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2743, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2746, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2749, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2752, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2755, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2758, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2761, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2764, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2767, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2770, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2773, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2776, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2785, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2794, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2803, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2812, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2821, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2830, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2839, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2848, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2857, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2866, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2875, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2884, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2893, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2902, 3, 25, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2920, 3, 28, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2938, 3, 31, 85,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2956, 3, 34, 91,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2974, 3, 37, 97,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2992, 3, 40, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3010, 3, 43, 109,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3028, 3, 46, 115,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3046, 3, 49, 121,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3064, 3, 52, 127,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3082, 3, 55, 133,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3100, 3, 58, 139,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3118, 3, 61, 145,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3136, 3, 64, 151,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3154, 3, 67, 157,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3172, 3, 73, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3202, 3, 79, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3232, 3, 85, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3262, 3, 91, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3292, 3, 97, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3322, 3, 103, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3352, 3, 109, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3382, 3, 115, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3412, 3, 121, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3442, 3, 127, 253,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3472, 3, 133, 263,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3502, 3, 139, 273,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3532, 3, 145, 283,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3562, 3, 151, 293,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3592, 3, 163, 303,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3637, 3, 173, 318,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3682, 3, 183, 333,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3727, 3, 193, 348,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3772, 3, 203, 363,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3817, 3, 213, 378,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3862, 3, 223, 393,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3907, 3, 233, 408,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3952, 3, 243, 423,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3997, 3, 253, 438,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4042, 3, 263, 453,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4087, 3, 273, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4132, 3, 283, 483,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4177, 3, 303, 498,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4240, 3, 318, 519,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4303, 3, 333, 540,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4366, 3, 348, 561,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4429, 3, 363, 582,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4492, 3, 378, 603,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4555, 3, 393, 624,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4618, 3, 408, 645,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4681, 3, 423, 666,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4744, 3, 438, 687,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4807, 3, 453, 708,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4870, 3, 468, 729,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4933, 3, 498, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5017, 3, 519, 778,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5101, 3, 540, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5185, 3, 561, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5269, 3, 582, 862,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5353, 3, 603, 890,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5437, 3, 624, 918,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5521, 3, 645, 946,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5605, 3, 666, 974,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5689, 3, 687,
                                                                       1002, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5773, 3, 708,
                                                                       1030, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5857, 3, 750,
                                                                       1058, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5965, 3, 778,
                                                                       1094, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6073, 3, 806,
                                                                       1130, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6181, 3, 834,
                                                                       1166, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6289, 3, 862,
                                                                       1202, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6397, 3, 890,
                                                                       1238, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6505, 3, 918,
                                                                       1274, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6613, 3, 946,
                                                                       1310, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6721, 3, 974,
                                                                       1346, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6829, 3, 1002,
                                                                       1382, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6937, 3, 1058,
                                                                       1418, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7072, 3, 1094,
                                                                       1463, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7207, 3, 1130,
                                                                       1508, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7342, 3, 1166,
                                                                       1553, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7477, 3, 1202,
                                                                       1598, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7612, 3, 1238,
                                                                       1643, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7747, 3, 1274,
                                                                       1688, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7882, 3, 1310,
                                                                       1733, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 8017, 3, 1346,
                                                                       1778, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8152, 3, 1418,
                                                                       1823, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8317, 3, 1463,
                                                                       1878, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8482, 3, 1508,
                                                                       1933, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8647, 3, 1553,
                                                                       1988, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8812, 3, 1598,
                                                                       2043, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8977, 3, 1643,
                                                                       2098, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9142, 3, 1688,
                                                                       2153, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 9307, 3, 1733,
                                                                       2208, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9472, 3, 1823,
                                                                       2263, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9670, 3, 1878,
                                                                       2329, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9868, 3, 1933,
                                                                       2395, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10066, 3, 1988,
                                                                       2461, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10264, 3, 2043,
                                                                       2527, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10462, 3, 2098,
                                                                       2593, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 10660, 3, 2153,
                                                                       2659, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10858, 3, 8, 9,
                                                                       2731, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10864, 3, 9, 10,
                                                                       2734, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10870, 3, 10, 11,
                                                                       2737, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10876, 3, 11, 12,
                                                                       2740, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10882, 3, 12, 13,
                                                                       2743, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10888, 3, 13, 14,
                                                                       2746, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10894, 3, 14, 15,
                                                                       2749, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10900, 3, 15, 16,
                                                                       2752, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10906, 3, 16, 17,
                                                                       2755, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10912, 3, 17, 18,
                                                                       2758, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10918, 3, 18, 19,
                                                                       2761, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10924, 3, 19, 20,
                                                                       2764, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10930, 3, 20, 21,
                                                                       2767, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10936, 3, 21, 22,
                                                                       2770, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10942, 3, 22, 23,
                                                                       2773, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10948, 0, 3,
                                                                       10858, 2731, 10864, 2776,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10966, 0, 3,
                                                                       10864, 2734, 10870, 2785,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 10984, 0, 3,
                                                                       10870, 2737, 10876, 2794,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11002, 0, 3,
                                                                       10876, 2740, 10882, 2803,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11020, 0, 3,
                                                                       10882, 2743, 10888, 2812,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11038, 0, 3,
                                                                       10888, 2746, 10894, 2821,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11056, 0, 3,
                                                                       10894, 2749, 10900, 2830,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11074, 0, 3,
                                                                       10900, 2752, 10906, 2839,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11092, 0, 3,
                                                                       10906, 2755, 10912, 2848,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11110, 0, 3,
                                                                       10912, 2758, 10918, 2857,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11128, 0, 3,
                                                                       10918, 2761, 10924, 2866,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11146, 0, 3,
                                                                       10924, 2764, 10930, 2875,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11164, 0, 3,
                                                                       10930, 2767, 10936, 2884,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11182, 0, 3,
                                                                       10936, 2770, 10942, 2893,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11200, 0, 3,
                                                                       10948, 2776, 10966, 73,
                                                                       79, 2938, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11236, 0, 3,
                                                                       10966, 2785, 10984, 79,
                                                                       85, 2956, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11272, 0, 3,
                                                                       10984, 2794, 11002, 85,
                                                                       91, 2974, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11308, 0, 3,
                                                                       11002, 2803, 11020, 91,
                                                                       97, 2992, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11344, 0, 3,
                                                                       11020, 2812, 11038, 97,
                                                                       103, 3010, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11380, 0, 3,
                                                                       11038, 2821, 11056, 103,
                                                                       109, 3028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11416, 0, 3,
                                                                       11056, 2830, 11074, 109,
                                                                       115, 3046, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11452, 0, 3,
                                                                       11074, 2839, 11092, 115,
                                                                       121, 3064, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11488, 0, 3,
                                                                       11092, 2848, 11110, 121,
                                                                       127, 3082, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11524, 0, 3,
                                                                       11110, 2857, 11128, 127,
                                                                       133, 3100, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11560, 0, 3,
                                                                       11128, 2866, 11146, 133,
                                                                       139, 3118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11596, 0, 3,
                                                                       11146, 2875, 11164, 139,
                                                                       145, 3136, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11632, 0, 3,
                                                                       11164, 2884, 11182, 145,
                                                                       151, 3154, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11668, 0, 3,
                                                                       11200, 2938, 11236, 163,
                                                                       173, 3232, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11728, 0, 3,
                                                                       11236, 2956, 11272, 173,
                                                                       183, 3262, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11788, 0, 3,
                                                                       11272, 2974, 11308, 183,
                                                                       193, 3292, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11848, 0, 3,
                                                                       11308, 2992, 11344, 193,
                                                                       203, 3322, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11908, 0, 3,
                                                                       11344, 3010, 11380, 203,
                                                                       213, 3352, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 11968, 0, 3,
                                                                       11380, 3028, 11416, 213,
                                                                       223, 3382, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12028, 0, 3,
                                                                       11416, 3046, 11452, 223,
                                                                       233, 3412, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12088, 0, 3,
                                                                       11452, 3064, 11488, 233,
                                                                       243, 3442, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12148, 0, 3,
                                                                       11488, 3082, 11524, 243,
                                                                       253, 3472, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12208, 0, 3,
                                                                       11524, 3100, 11560, 253,
                                                                       263, 3502, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12268, 0, 3,
                                                                       11560, 3118, 11596, 263,
                                                                       273, 3532, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12328, 0, 3,
                                                                       11596, 3136, 11632, 273,
                                                                       283, 3562, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12388, 0, 3,
                                                                       11668, 3232, 11728, 303,
                                                                       318, 3682, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12478, 0, 3,
                                                                       11728, 3262, 11788, 318,
                                                                       333, 3727, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12568, 0, 3,
                                                                       11788, 3292, 11848, 333,
                                                                       348, 3772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12658, 0, 3,
                                                                       11848, 3322, 11908, 348,
                                                                       363, 3817, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12748, 0, 3,
                                                                       11908, 3352, 11968, 363,
                                                                       378, 3862, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12838, 0, 3,
                                                                       11968, 3382, 12028, 378,
                                                                       393, 3907, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12928, 0, 3,
                                                                       12028, 3412, 12088, 393,
                                                                       408, 3952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13018, 0, 3,
                                                                       12088, 3442, 12148, 408,
                                                                       423, 3997, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13108, 0, 3,
                                                                       12148, 3472, 12208, 423,
                                                                       438, 4042, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13198, 0, 3,
                                                                       12208, 3502, 12268, 438,
                                                                       453, 4087, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13288, 0, 3,
                                                                       12268, 3532, 12328, 453,
                                                                       468, 4132, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13378, 0, 3,
                                                                       12388, 3682, 12478, 498,
                                                                       519, 4303, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13504, 0, 3,
                                                                       12478, 3727, 12568, 519,
                                                                       540, 4366, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13630, 0, 3,
                                                                       12568, 3772, 12658, 540,
                                                                       561, 4429, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13756, 0, 3,
                                                                       12658, 3817, 12748, 561,
                                                                       582, 4492, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13882, 0, 3,
                                                                       12748, 3862, 12838, 582,
                                                                       603, 4555, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14008, 0, 3,
                                                                       12838, 3907, 12928, 603,
                                                                       624, 4618, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14134, 0, 3,
                                                                       12928, 3952, 13018, 624,
                                                                       645, 4681, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14260, 0, 3,
                                                                       13018, 3997, 13108, 645,
                                                                       666, 4744, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14386, 0, 3,
                                                                       13108, 4042, 13198, 666,
                                                                       687, 4807, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14512, 0, 3,
                                                                       13198, 4087, 13288, 687,
                                                                       708, 4870, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14638, 0, 3,
                                                                       13378, 4303, 13504, 750,
                                                                       778, 5101, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14806, 0, 3,
                                                                       13504, 4366, 13630, 778,
                                                                       806, 5185, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14974, 0, 3,
                                                                       13630, 4429, 13756, 806,
                                                                       834, 5269, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15142, 0, 3,
                                                                       13756, 4492, 13882, 834,
                                                                       862, 5353, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15310, 0, 3,
                                                                       13882, 4555, 14008, 862,
                                                                       890, 5437, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15478, 0, 3,
                                                                       14008, 4618, 14134, 890,
                                                                       918, 5521, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15646, 0, 3,
                                                                       14134, 4681, 14260, 918,
                                                                       946, 5605, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15814, 0, 3,
                                                                       14260, 4744, 14386, 946,
                                                                       974, 5689, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15982, 0, 3,
                                                                       14386, 4807, 14512, 974,
                                                                       1002, 5773, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16150, 0, 3,
                                                                       14638, 5101, 14806, 1058,
                                                                       1094, 6073, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16366, 0, 3,
                                                                       14806, 5185, 14974, 1094,
                                                                       1130, 6181, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16582, 0, 3,
                                                                       14974, 5269, 15142, 1130,
                                                                       1166, 6289, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16798, 0, 3,
                                                                       15142, 5353, 15310, 1166,
                                                                       1202, 6397, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17014, 0, 3,
                                                                       15310, 5437, 15478, 1202,
                                                                       1238, 6505, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17230, 0, 3,
                                                                       15478, 5521, 15646, 1238,
                                                                       1274, 6613, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17446, 0, 3,
                                                                       15646, 5605, 15814, 1274,
                                                                       1310, 6721, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17662, 0, 3,
                                                                       15814, 5689, 15982, 1310,
                                                                       1346, 6829, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17878, 0, 3,
                                                                       16150, 6073, 16366, 1418,
                                                                       1463, 7207, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18148, 0, 3,
                                                                       16366, 6181, 16582, 1463,
                                                                       1508, 7342, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18418, 0, 3,
                                                                       16582, 6289, 16798, 1508,
                                                                       1553, 7477, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18688, 0, 3,
                                                                       16798, 6397, 17014, 1553,
                                                                       1598, 7612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18958, 0, 3,
                                                                       17014, 6505, 17230, 1598,
                                                                       1643, 7747, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 19228, 0, 3,
                                                                       17230, 6613, 17446, 1643,
                                                                       1688, 7882, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 19498, 0, 3,
                                                                       17446, 6721, 17662, 1688,
                                                                       1733, 8017, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19768, 0, 3,
                                                                       17878, 7207, 18148, 1823,
                                                                       1878, 8482, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20098, 0, 3,
                                                                       18148, 7342, 18418, 1878,
                                                                       1933, 8647, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20428, 0, 3,
                                                                       18418, 7477, 18688, 1933,
                                                                       1988, 8812, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20758, 0, 3,
                                                                       18688, 7612, 18958, 1988,
                                                                       2043, 8977, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 21088, 0, 3,
                                                                       18958, 7747, 19228, 2043,
                                                                       2098, 9142, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 21418, 0, 3,
                                                                       19228, 7882, 19498, 2098,
                                                                       2153, 9307, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 21748, 0, 3,
                                                                       19768, 8482, 20098, 2263,
                                                                       2329, 9868, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 22144, 0, 3,
                                                                       20098, 8647, 20428, 2329,
                                                                       2395, 10066, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 22540, 0, 3,
                                                                       20428, 8812, 20758, 2395,
                                                                       2461, 10264, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 22936, 0, 3,
                                                                       20758, 8977, 21088, 2461,
                                                                       2527, 10462, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 23332, 0, 3,
                                                                       21088, 9142, 21418, 2527,
                                                                       2593, 10660, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23728, 3, 2725,
                                                                       2728, 10858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23738, 3, 2728,
                                                                       2731, 10864, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23748, 3, 2731,
                                                                       2734, 10870, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23758, 3, 2734,
                                                                       2737, 10876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23768, 3, 2737,
                                                                       2740, 10882, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23778, 3, 2740,
                                                                       2743, 10888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23788, 3, 2743,
                                                                       2746, 10894, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23798, 3, 2746,
                                                                       2749, 10900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23808, 3, 2749,
                                                                       2752, 10906, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23818, 3, 2752,
                                                                       2755, 10912, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23828, 3, 2755,
                                                                       2758, 10918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23838, 3, 2758,
                                                                       2761, 10924, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23848, 3, 2761,
                                                                       2764, 10930, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23858, 3, 2764,
                                                                       2767, 10936, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23868, 3, 2767,
                                                                       2770, 10942, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23878, 0, 3,
                                                                       23728, 10858, 23738,
                                                                       10948, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23908, 0, 3,
                                                                       23738, 10864, 23748,
                                                                       10966, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23938, 0, 3,
                                                                       23748, 10870, 23758,
                                                                       10984, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23968, 0, 3,
                                                                       23758, 10876, 23768,
                                                                       11002, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23998, 0, 3,
                                                                       23768, 10882, 23778,
                                                                       11020, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24028, 0, 3,
                                                                       23778, 10888, 23788,
                                                                       11038, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24058, 0, 3,
                                                                       23788, 10894, 23798,
                                                                       11056, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24088, 0, 3,
                                                                       23798, 10900, 23808,
                                                                       11074, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24118, 0, 3,
                                                                       23808, 10906, 23818,
                                                                       11092, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24148, 0, 3,
                                                                       23818, 10912, 23828,
                                                                       11110, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24178, 0, 3,
                                                                       23828, 10918, 23838,
                                                                       11128, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24208, 0, 3,
                                                                       23838, 10924, 23848,
                                                                       11146, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24238, 0, 3,
                                                                       23848, 10930, 23858,
                                                                       11164, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24268, 0, 3,
                                                                       23858, 10936, 23868,
                                                                       11182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24298, 0, 3,
                                                                       23878, 10948, 23908, 2902,
                                                                       2920, 11200, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24358, 0, 3,
                                                                       23908, 10966, 23938, 2920,
                                                                       2938, 11236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24418, 0, 3,
                                                                       23938, 10984, 23968, 2938,
                                                                       2956, 11272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24478, 0, 3,
                                                                       23968, 11002, 23998, 2956,
                                                                       2974, 11308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24538, 0, 3,
                                                                       23998, 11020, 24028, 2974,
                                                                       2992, 11344, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24598, 0, 3,
                                                                       24028, 11038, 24058, 2992,
                                                                       3010, 11380, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24658, 0, 3,
                                                                       24058, 11056, 24088, 3010,
                                                                       3028, 11416, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24718, 0, 3,
                                                                       24088, 11074, 24118, 3028,
                                                                       3046, 11452, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24778, 0, 3,
                                                                       24118, 11092, 24148, 3046,
                                                                       3064, 11488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24838, 0, 3,
                                                                       24148, 11110, 24178, 3064,
                                                                       3082, 11524, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24898, 0, 3,
                                                                       24178, 11128, 24208, 3082,
                                                                       3100, 11560, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24958, 0, 3,
                                                                       24208, 11146, 24238, 3100,
                                                                       3118, 11596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25018, 0, 3,
                                                                       24238, 11164, 24268, 3118,
                                                                       3136, 11632, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25078, 0, 3,
                                                                       24298, 11200, 24358, 3172,
                                                                       3202, 11668, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25178, 0, 3,
                                                                       24358, 11236, 24418, 3202,
                                                                       3232, 11728, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25278, 0, 3,
                                                                       24418, 11272, 24478, 3232,
                                                                       3262, 11788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25378, 0, 3,
                                                                       24478, 11308, 24538, 3262,
                                                                       3292, 11848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25478, 0, 3,
                                                                       24538, 11344, 24598, 3292,
                                                                       3322, 11908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25578, 0, 3,
                                                                       24598, 11380, 24658, 3322,
                                                                       3352, 11968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25678, 0, 3,
                                                                       24658, 11416, 24718, 3352,
                                                                       3382, 12028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25778, 0, 3,
                                                                       24718, 11452, 24778, 3382,
                                                                       3412, 12088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25878, 0, 3,
                                                                       24778, 11488, 24838, 3412,
                                                                       3442, 12148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25978, 0, 3,
                                                                       24838, 11524, 24898, 3442,
                                                                       3472, 12208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26078, 0, 3,
                                                                       24898, 11560, 24958, 3472,
                                                                       3502, 12268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26178, 0, 3,
                                                                       24958, 11596, 25018, 3502,
                                                                       3532, 12328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26278, 0, 3,
                                                                       25078, 11668, 25178, 3592,
                                                                       3637, 12388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26428, 0, 3,
                                                                       25178, 11728, 25278, 3637,
                                                                       3682, 12478, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26578, 0, 3,
                                                                       25278, 11788, 25378, 3682,
                                                                       3727, 12568, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26728, 0, 3,
                                                                       25378, 11848, 25478, 3727,
                                                                       3772, 12658, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26878, 0, 3,
                                                                       25478, 11908, 25578, 3772,
                                                                       3817, 12748, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27028, 0, 3,
                                                                       25578, 11968, 25678, 3817,
                                                                       3862, 12838, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27178, 0, 3,
                                                                       25678, 12028, 25778, 3862,
                                                                       3907, 12928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27328, 0, 3,
                                                                       25778, 12088, 25878, 3907,
                                                                       3952, 13018, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27478, 0, 3,
                                                                       25878, 12148, 25978, 3952,
                                                                       3997, 13108, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27628, 0, 3,
                                                                       25978, 12208, 26078, 3997,
                                                                       4042, 13198, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27778, 0, 3,
                                                                       26078, 12268, 26178, 4042,
                                                                       4087, 13288, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27928, 0, 3,
                                                                       26278, 12388, 26428, 4177,
                                                                       4240, 13378, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28138, 0, 3,
                                                                       26428, 12478, 26578, 4240,
                                                                       4303, 13504, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28348, 0, 3,
                                                                       26578, 12568, 26728, 4303,
                                                                       4366, 13630, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28558, 0, 3,
                                                                       26728, 12658, 26878, 4366,
                                                                       4429, 13756, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28768, 0, 3,
                                                                       26878, 12748, 27028, 4429,
                                                                       4492, 13882, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28978, 0, 3,
                                                                       27028, 12838, 27178, 4492,
                                                                       4555, 14008, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29188, 0, 3,
                                                                       27178, 12928, 27328, 4555,
                                                                       4618, 14134, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29398, 0, 3,
                                                                       27328, 13018, 27478, 4618,
                                                                       4681, 14260, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29608, 0, 3,
                                                                       27478, 13108, 27628, 4681,
                                                                       4744, 14386, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29818, 0, 3,
                                                                       27628, 13198, 27778, 4744,
                                                                       4807, 14512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30028, 0, 3,
                                                                       27928, 13378, 28138, 4933,
                                                                       5017, 14638, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30308, 0, 3,
                                                                       28138, 13504, 28348, 5017,
                                                                       5101, 14806, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30588, 0, 3,
                                                                       28348, 13630, 28558, 5101,
                                                                       5185, 14974, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30868, 0, 3,
                                                                       28558, 13756, 28768, 5185,
                                                                       5269, 15142, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31148, 0, 3,
                                                                       28768, 13882, 28978, 5269,
                                                                       5353, 15310, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31428, 0, 3,
                                                                       28978, 14008, 29188, 5353,
                                                                       5437, 15478, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31708, 0, 3,
                                                                       29188, 14134, 29398, 5437,
                                                                       5521, 15646, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31988, 0, 3,
                                                                       29398, 14260, 29608, 5521,
                                                                       5605, 15814, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 32268, 0, 3,
                                                                       29608, 14386, 29818, 5605,
                                                                       5689, 15982, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32548, 0, 3,
                                                                       30028, 14638, 30308, 5857,
                                                                       5965, 16150, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32908, 0, 3,
                                                                       30308, 14806, 30588, 5965,
                                                                       6073, 16366, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33268, 0, 3,
                                                                       30588, 14974, 30868, 6073,
                                                                       6181, 16582, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33628, 0, 3,
                                                                       30868, 15142, 31148, 6181,
                                                                       6289, 16798, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33988, 0, 3,
                                                                       31148, 15310, 31428, 6289,
                                                                       6397, 17014, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 34348, 0, 3,
                                                                       31428, 15478, 31708, 6397,
                                                                       6505, 17230, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 34708, 0, 3,
                                                                       31708, 15646, 31988, 6505,
                                                                       6613, 17446, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 35068, 0, 3,
                                                                       31988, 15814, 32268, 6613,
                                                                       6721, 17662, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35428, 0, 3,
                                                                       32548, 16150, 32908, 6937,
                                                                       7072, 17878, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35878, 0, 3,
                                                                       32908, 16366, 33268, 7072,
                                                                       7207, 18148, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 36328, 0, 3,
                                                                       33268, 16582, 33628, 7207,
                                                                       7342, 18418, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 36778, 0, 3,
                                                                       33628, 16798, 33988, 7342,
                                                                       7477, 18688, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 37228, 0, 3,
                                                                       33988, 17014, 34348, 7477,
                                                                       7612, 18958, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 37678, 0, 3,
                                                                       34348, 17230, 34708, 7612,
                                                                       7747, 19228, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 38128, 0, 3,
                                                                       34708, 17446, 35068, 7747,
                                                                       7882, 19498, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 38578, 0, 3,
                                                                       35428, 17878, 35878, 8152,
                                                                       8317, 19768, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 39128, 0, 3,
                                                                       35878, 18148, 36328, 8317,
                                                                       8482, 20098, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 39678, 0, 3,
                                                                       36328, 18418, 36778, 8482,
                                                                       8647, 20428, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 40228, 0, 3,
                                                                       36778, 18688, 37228, 8647,
                                                                       8812, 20758, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 40778, 0, 3,
                                                                       37228, 18958, 37678, 8812,
                                                                       8977, 21088, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 41328, 0, 3,
                                                                       37678, 19228, 38128, 8977,
                                                                       9142, 21418, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 41878, 0, 3,
                                                                       38578, 19768, 39128, 9472,
                                                                       9670, 21748, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 42538, 0, 3,
                                                                       39128, 20098, 39678, 9670,
                                                                       9868, 22144, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 43198, 0, 3,
                                                                       39678, 20428, 40228, 9868,
                                                                       10066, 22540, ncols,
                                                                       gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 43858, 0, 3,
                                                                       40228, 20758, 40778,
                                                                       10066, 10264, 22936,
                                                                       ncols, gamma, p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 44518, 0, 3,
                                                                       40778, 21088, 41328,
                                                                       10264, 10462, 23332,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45178, 3, 10858,
                                                                       10864, 23748, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45193, 3, 10864,
                                                                       10870, 23758, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45208, 3, 10870,
                                                                       10876, 23768, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45223, 3, 10876,
                                                                       10882, 23778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45238, 3, 10882,
                                                                       10888, 23788, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45253, 3, 10888,
                                                                       10894, 23798, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45268, 3, 10894,
                                                                       10900, 23808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45283, 3, 10900,
                                                                       10906, 23818, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45298, 3, 10906,
                                                                       10912, 23828, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45313, 3, 10912,
                                                                       10918, 23838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45328, 3, 10918,
                                                                       10924, 23848, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45343, 3, 10924,
                                                                       10930, 23858, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 45358, 3, 10930,
                                                                       10936, 23868, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45373, 0, 3,
                                                                       45178, 23748, 45193,
                                                                       10948, 10966, 23938,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45418, 0, 3,
                                                                       45193, 23758, 45208,
                                                                       10966, 10984, 23968,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45463, 0, 3,
                                                                       45208, 23768, 45223,
                                                                       10984, 11002, 23998,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45508, 0, 3,
                                                                       45223, 23778, 45238,
                                                                       11002, 11020, 24028,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45553, 0, 3,
                                                                       45238, 23788, 45253,
                                                                       11020, 11038, 24058,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45598, 0, 3,
                                                                       45253, 23798, 45268,
                                                                       11038, 11056, 24088,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45643, 0, 3,
                                                                       45268, 23808, 45283,
                                                                       11056, 11074, 24118,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45688, 0, 3,
                                                                       45283, 23818, 45298,
                                                                       11074, 11092, 24148,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45733, 0, 3,
                                                                       45298, 23828, 45313,
                                                                       11092, 11110, 24178,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45778, 0, 3,
                                                                       45313, 23838, 45328,
                                                                       11110, 11128, 24208,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45823, 0, 3,
                                                                       45328, 23848, 45343,
                                                                       11128, 11146, 24238,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45868, 0, 3,
                                                                       45343, 23858, 45358,
                                                                       11146, 11164, 24268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45913, 0, 3,
                                                                       45373, 23938, 45418,
                                                                       11200, 11236, 24418,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46003, 0, 3,
                                                                       45418, 23968, 45463,
                                                                       11236, 11272, 24478,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46093, 0, 3,
                                                                       45463, 23998, 45508,
                                                                       11272, 11308, 24538,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46183, 0, 3,
                                                                       45508, 24028, 45553,
                                                                       11308, 11344, 24598,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46273, 0, 3,
                                                                       45553, 24058, 45598,
                                                                       11344, 11380, 24658,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46363, 0, 3,
                                                                       45598, 24088, 45643,
                                                                       11380, 11416, 24718,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46453, 0, 3,
                                                                       45643, 24118, 45688,
                                                                       11416, 11452, 24778,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46543, 0, 3,
                                                                       45688, 24148, 45733,
                                                                       11452, 11488, 24838,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46633, 0, 3,
                                                                       45733, 24178, 45778,
                                                                       11488, 11524, 24898,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46723, 0, 3,
                                                                       45778, 24208, 45823,
                                                                       11524, 11560, 24958,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46813, 0, 3,
                                                                       45823, 24238, 45868,
                                                                       11560, 11596, 25018,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46903, 0, 3,
                                                                       45913, 24418, 46003,
                                                                       11668, 11728, 25278,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47053, 0, 3,
                                                                       46003, 24478, 46093,
                                                                       11728, 11788, 25378,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47203, 0, 3,
                                                                       46093, 24538, 46183,
                                                                       11788, 11848, 25478,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47353, 0, 3,
                                                                       46183, 24598, 46273,
                                                                       11848, 11908, 25578,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47503, 0, 3,
                                                                       46273, 24658, 46363,
                                                                       11908, 11968, 25678,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47653, 0, 3,
                                                                       46363, 24718, 46453,
                                                                       11968, 12028, 25778,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47803, 0, 3,
                                                                       46453, 24778, 46543,
                                                                       12028, 12088, 25878,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47953, 0, 3,
                                                                       46543, 24838, 46633,
                                                                       12088, 12148, 25978,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 48103, 0, 3,
                                                                       46633, 24898, 46723,
                                                                       12148, 12208, 26078,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 48253, 0, 3,
                                                                       46723, 24958, 46813,
                                                                       12208, 12268, 26178,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48403, 0, 3,
                                                                       46903, 25278, 47053,
                                                                       12388, 12478, 26578,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48628, 0, 3,
                                                                       47053, 25378, 47203,
                                                                       12478, 12568, 26728,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48853, 0, 3,
                                                                       47203, 25478, 47353,
                                                                       12568, 12658, 26878,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49078, 0, 3,
                                                                       47353, 25578, 47503,
                                                                       12658, 12748, 27028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49303, 0, 3,
                                                                       47503, 25678, 47653,
                                                                       12748, 12838, 27178,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49528, 0, 3,
                                                                       47653, 25778, 47803,
                                                                       12838, 12928, 27328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49753, 0, 3,
                                                                       47803, 25878, 47953,
                                                                       12928, 13018, 27478,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49978, 0, 3,
                                                                       47953, 25978, 48103,
                                                                       13018, 13108, 27628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 50203, 0, 3,
                                                                       48103, 26078, 48253,
                                                                       13108, 13198, 27778,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50428, 0, 3,
                                                                       48403, 26578, 48628,
                                                                       13378, 13504, 28348,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50743, 0, 3,
                                                                       48628, 26728, 48853,
                                                                       13504, 13630, 28558,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51058, 0, 3,
                                                                       48853, 26878, 49078,
                                                                       13630, 13756, 28768,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51373, 0, 3,
                                                                       49078, 27028, 49303,
                                                                       13756, 13882, 28978,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51688, 0, 3,
                                                                       49303, 27178, 49528,
                                                                       13882, 14008, 29188,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 52003, 0, 3,
                                                                       49528, 27328, 49753,
                                                                       14008, 14134, 29398,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 52318, 0, 3,
                                                                       49753, 27478, 49978,
                                                                       14134, 14260, 29608,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 52633, 0, 3,
                                                                       49978, 27628, 50203,
                                                                       14260, 14386, 29818,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 52948, 0, 3,
                                                                       50428, 28348, 50743,
                                                                       14638, 14806, 30588,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 53368, 0, 3,
                                                                       50743, 28558, 51058,
                                                                       14806, 14974, 30868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 53788, 0, 3,
                                                                       51058, 28768, 51373,
                                                                       14974, 15142, 31148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 54208, 0, 3,
                                                                       51373, 28978, 51688,
                                                                       15142, 15310, 31428,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 54628, 0, 3,
                                                                       51688, 29188, 52003,
                                                                       15310, 15478, 31708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 55048, 0, 3,
                                                                       52003, 29398, 52318,
                                                                       15478, 15646, 31988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 55468, 0, 3,
                                                                       52318, 29608, 52633,
                                                                       15646, 15814, 32268,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 55888, 0, 3,
                                                                       52948, 30588, 53368,
                                                                       16150, 16366, 33268,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 56428, 0, 3,
                                                                       53368, 30868, 53788,
                                                                       16366, 16582, 33628,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 56968, 0, 3,
                                                                       53788, 31148, 54208,
                                                                       16582, 16798, 33988,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 57508, 0, 3,
                                                                       54208, 31428, 54628,
                                                                       16798, 17014, 34348,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 58048, 0, 3,
                                                                       54628, 31708, 55048,
                                                                       17014, 17230, 34708,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 58588, 0, 3,
                                                                       55048, 31988, 55468,
                                                                       17230, 17446, 35068,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 59128, 0, 3,
                                                                       55888, 33268, 56428,
                                                                       17878, 18148, 36328,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 59803, 0, 3,
                                                                       56428, 33628, 56968,
                                                                       18148, 18418, 36778,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 60478, 0, 3,
                                                                       56968, 33988, 57508,
                                                                       18418, 18688, 37228,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 61153, 0, 3,
                                                                       57508, 34348, 58048,
                                                                       18688, 18958, 37678,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 61828, 0, 3,
                                                                       58048, 34708, 58588,
                                                                       18958, 19228, 38128,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 62503, 0, 3,
                                                                       59128, 36328, 59803,
                                                                       19768, 20098, 39678,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 63328, 0, 3,
                                                                       59803, 36778, 60478,
                                                                       20098, 20428, 40228,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 64153, 0, 3,
                                                                       60478, 37228, 61153,
                                                                       20428, 20758, 40778,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 64978, 0, 3,
                                                                       61153, 37678, 61828,
                                                                       20758, 21088, 41328,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 65803, 0, 3,
                                                                       62503, 39678, 63328,
                                                                       21748, 22144, 43198,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 66793, 0, 3,
                                                                       63328, 40228, 64153,
                                                                       22144, 22540, 43858,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 67783, 0, 3,
                                                                       64153, 40778, 64978,
                                                                       22540, 22936, 44518,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68773, 3, 23728,
                                                                       23738, 45178, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68794, 3, 23738,
                                                                       23748, 45193, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68815, 3, 23748,
                                                                       23758, 45208, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68836, 3, 23758,
                                                                       23768, 45223, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68857, 3, 23768,
                                                                       23778, 45238, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68878, 3, 23778,
                                                                       23788, 45253, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68899, 3, 23788,
                                                                       23798, 45268, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68920, 3, 23798,
                                                                       23808, 45283, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68941, 3, 23808,
                                                                       23818, 45298, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68962, 3, 23818,
                                                                       23828, 45313, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 68983, 3, 23828,
                                                                       23838, 45328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 69004, 3, 23838,
                                                                       23848, 45343, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 69025, 3, 23848,
                                                                       23858, 45358, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69046, 0, 3,
                                                                       68773, 45178, 68794,
                                                                       23878, 23908, 45373,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69109, 0, 3,
                                                                       68794, 45193, 68815,
                                                                       23908, 23938, 45418,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69172, 0, 3,
                                                                       68815, 45208, 68836,
                                                                       23938, 23968, 45463,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69235, 0, 3,
                                                                       68836, 45223, 68857,
                                                                       23968, 23998, 45508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69298, 0, 3,
                                                                       68857, 45238, 68878,
                                                                       23998, 24028, 45553,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69361, 0, 3,
                                                                       68878, 45253, 68899,
                                                                       24028, 24058, 45598,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69424, 0, 3,
                                                                       68899, 45268, 68920,
                                                                       24058, 24088, 45643,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69487, 0, 3,
                                                                       68920, 45283, 68941,
                                                                       24088, 24118, 45688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69550, 0, 3,
                                                                       68941, 45298, 68962,
                                                                       24118, 24148, 45733,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69613, 0, 3,
                                                                       68962, 45313, 68983,
                                                                       24148, 24178, 45778,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69676, 0, 3,
                                                                       68983, 45328, 69004,
                                                                       24178, 24208, 45823,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 69739, 0, 3,
                                                                       69004, 45343, 69025,
                                                                       24208, 24238, 45868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 69802, 0, 3,
                                                                       69046, 45373, 69109,
                                                                       24298, 24358, 45913,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 69928, 0, 3,
                                                                       69109, 45418, 69172,
                                                                       24358, 24418, 46003,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70054, 0, 3,
                                                                       69172, 45463, 69235,
                                                                       24418, 24478, 46093,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70180, 0, 3,
                                                                       69235, 45508, 69298,
                                                                       24478, 24538, 46183,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70306, 0, 3,
                                                                       69298, 45553, 69361,
                                                                       24538, 24598, 46273,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70432, 0, 3,
                                                                       69361, 45598, 69424,
                                                                       24598, 24658, 46363,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70558, 0, 3,
                                                                       69424, 45643, 69487,
                                                                       24658, 24718, 46453,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70684, 0, 3,
                                                                       69487, 45688, 69550,
                                                                       24718, 24778, 46543,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70810, 0, 3,
                                                                       69550, 45733, 69613,
                                                                       24778, 24838, 46633,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 70936, 0, 3,
                                                                       69613, 45778, 69676,
                                                                       24838, 24898, 46723,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 71062, 0, 3,
                                                                       69676, 45823, 69739,
                                                                       24898, 24958, 46813,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 71188, 0, 3,
                                                                       69802, 45913, 69928,
                                                                       25078, 25178, 46903,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 71398, 0, 3,
                                                                       69928, 46003, 70054,
                                                                       25178, 25278, 47053,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 71608, 0, 3,
                                                                       70054, 46093, 70180,
                                                                       25278, 25378, 47203,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 71818, 0, 3,
                                                                       70180, 46183, 70306,
                                                                       25378, 25478, 47353,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 72028, 0, 3,
                                                                       70306, 46273, 70432,
                                                                       25478, 25578, 47503,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 72238, 0, 3,
                                                                       70432, 46363, 70558,
                                                                       25578, 25678, 47653,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 72448, 0, 3,
                                                                       70558, 46453, 70684,
                                                                       25678, 25778, 47803,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 72658, 0, 3,
                                                                       70684, 46543, 70810,
                                                                       25778, 25878, 47953,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 72868, 0, 3,
                                                                       70810, 46633, 70936,
                                                                       25878, 25978, 48103,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 73078, 0, 3,
                                                                       70936, 46723, 71062,
                                                                       25978, 26078, 48253,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 73288, 0, 3,
                                                                       71188, 46903, 71398,
                                                                       26278, 26428, 48403,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 73603, 0, 3,
                                                                       71398, 47053, 71608,
                                                                       26428, 26578, 48628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 73918, 0, 3,
                                                                       71608, 47203, 71818,
                                                                       26578, 26728, 48853,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 74233, 0, 3,
                                                                       71818, 47353, 72028,
                                                                       26728, 26878, 49078,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 74548, 0, 3,
                                                                       72028, 47503, 72238,
                                                                       26878, 27028, 49303,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 74863, 0, 3,
                                                                       72238, 47653, 72448,
                                                                       27028, 27178, 49528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 75178, 0, 3,
                                                                       72448, 47803, 72658,
                                                                       27178, 27328, 49753,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 75493, 0, 3,
                                                                       72658, 47953, 72868,
                                                                       27328, 27478, 49978,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 75808, 0, 3,
                                                                       72868, 48103, 73078,
                                                                       27478, 27628, 50203,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 76123, 0, 3,
                                                                       73288, 48403, 73603,
                                                                       27928, 28138, 50428,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 76564, 0, 3,
                                                                       73603, 48628, 73918,
                                                                       28138, 28348, 50743,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 77005, 0, 3,
                                                                       73918, 48853, 74233,
                                                                       28348, 28558, 51058,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 77446, 0, 3,
                                                                       74233, 49078, 74548,
                                                                       28558, 28768, 51373,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 77887, 0, 3,
                                                                       74548, 49303, 74863,
                                                                       28768, 28978, 51688,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 78328, 0, 3,
                                                                       74863, 49528, 75178,
                                                                       28978, 29188, 52003,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 78769, 0, 3,
                                                                       75178, 49753, 75493,
                                                                       29188, 29398, 52318,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 79210, 0, 3,
                                                                       75493, 49978, 75808,
                                                                       29398, 29608, 52633,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 79651, 0, 3,
                                                                       76123, 50428, 76564,
                                                                       30028, 30308, 52948,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 80239, 0, 3,
                                                                       76564, 50743, 77005,
                                                                       30308, 30588, 53368,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 80827, 0, 3,
                                                                       77005, 51058, 77446,
                                                                       30588, 30868, 53788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 81415, 0, 3,
                                                                       77446, 51373, 77887,
                                                                       30868, 31148, 54208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 82003, 0, 3,
                                                                       77887, 51688, 78328,
                                                                       31148, 31428, 54628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 82591, 0, 3,
                                                                       78328, 52003, 78769,
                                                                       31428, 31708, 55048,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 83179, 0, 3,
                                                                       78769, 52318, 79210,
                                                                       31708, 31988, 55468,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 83767, 0, 3,
                                                                       79651, 52948, 80239,
                                                                       32548, 32908, 55888,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 84523, 0, 3,
                                                                       80239, 53368, 80827,
                                                                       32908, 33268, 56428,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 85279, 0, 3,
                                                                       80827, 53788, 81415,
                                                                       33268, 33628, 56968,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 86035, 0, 3,
                                                                       81415, 54208, 82003,
                                                                       33628, 33988, 57508,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 86791, 0, 3,
                                                                       82003, 54628, 82591,
                                                                       33988, 34348, 58048,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 87547, 0, 3,
                                                                       82591, 55048, 83179,
                                                                       34348, 34708, 58588,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 88303, 0, 3,
                                                                       83767, 55888, 84523,
                                                                       35428, 35878, 59128,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 89248, 0, 3,
                                                                       84523, 56428, 85279,
                                                                       35878, 36328, 59803,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 90193, 0, 3,
                                                                       85279, 56968, 86035,
                                                                       36328, 36778, 60478,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 91138, 0, 3,
                                                                       86035, 57508, 86791,
                                                                       36778, 37228, 61153,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 92083, 0, 3,
                                                                       86791, 58048, 87547,
                                                                       37228, 37678, 61828,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 93028, 0, 3,
                                                                       88303, 59128, 89248,
                                                                       38578, 39128, 62503,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 94183, 0, 3,
                                                                       89248, 59803, 90193,
                                                                       39128, 39678, 63328,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 95338, 0, 3,
                                                                       90193, 60478, 91138,
                                                                       39678, 40228, 64153,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 96493, 0, 3,
                                                                       91138, 61153, 92083,
                                                                       40228, 40778, 64978,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 97648, 0, 3,
                                                                       93028, 62503, 94183,
                                                                       41878, 42538, 65803,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 99034, 0, 3,
                                                                       94183, 63328, 95338,
                                                                       42538, 43198, 66793,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 100420, 0, 3,
                                                                       95338, 64153, 96493,
                                                                       43198, 43858, 67783,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101806, 3, 45178,
                                                                       45193, 68815, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101834, 3, 45193,
                                                                       45208, 68836, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101862, 3, 45208,
                                                                       45223, 68857, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101890, 3, 45223,
                                                                       45238, 68878, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101918, 3, 45238,
                                                                       45253, 68899, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101946, 3, 45253,
                                                                       45268, 68920, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 101974, 3, 45268,
                                                                       45283, 68941, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102002, 3, 45283,
                                                                       45298, 68962, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102030, 3, 45298,
                                                                       45313, 68983, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102058, 3, 45313,
                                                                       45328, 69004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 102086, 3, 45328,
                                                                       45343, 69025, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102114, 0, 3,
                                                                       101806, 68815, 101834,
                                                                       45373, 45418, 69172,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102198, 0, 3,
                                                                       101834, 68836, 101862,
                                                                       45418, 45463, 69235,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102282, 0, 3,
                                                                       101862, 68857, 101890,
                                                                       45463, 45508, 69298,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102366, 0, 3,
                                                                       101890, 68878, 101918,
                                                                       45508, 45553, 69361,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102450, 0, 3,
                                                                       101918, 68899, 101946,
                                                                       45553, 45598, 69424,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102534, 0, 3,
                                                                       101946, 68920, 101974,
                                                                       45598, 45643, 69487,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102618, 0, 3,
                                                                       101974, 68941, 102002,
                                                                       45643, 45688, 69550,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102702, 0, 3,
                                                                       102002, 68962, 102030,
                                                                       45688, 45733, 69613,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102786, 0, 3,
                                                                       102030, 68983, 102058,
                                                                       45733, 45778, 69676,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 102870, 0, 3,
                                                                       102058, 69004, 102086,
                                                                       45778, 45823, 69739,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 102954, 0, 3,
                                                                       102114, 69172, 102198,
                                                                       45913, 46003, 70054,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103122, 0, 3,
                                                                       102198, 69235, 102282,
                                                                       46003, 46093, 70180,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103290, 0, 3,
                                                                       102282, 69298, 102366,
                                                                       46093, 46183, 70306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103458, 0, 3,
                                                                       102366, 69361, 102450,
                                                                       46183, 46273, 70432,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103626, 0, 3,
                                                                       102450, 69424, 102534,
                                                                       46273, 46363, 70558,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103794, 0, 3,
                                                                       102534, 69487, 102618,
                                                                       46363, 46453, 70684,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 103962, 0, 3,
                                                                       102618, 69550, 102702,
                                                                       46453, 46543, 70810,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 104130, 0, 3,
                                                                       102702, 69613, 102786,
                                                                       46543, 46633, 70936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 104298, 0, 3,
                                                                       102786, 69676, 102870,
                                                                       46633, 46723, 71062,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 104466, 0, 3,
                                                                       102954, 70054, 103122,
                                                                       46903, 47053, 71608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 104746, 0, 3,
                                                                       103122, 70180, 103290,
                                                                       47053, 47203, 71818,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 105026, 0, 3,
                                                                       103290, 70306, 103458,
                                                                       47203, 47353, 72028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 105306, 0, 3,
                                                                       103458, 70432, 103626,
                                                                       47353, 47503, 72238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 105586, 0, 3,
                                                                       103626, 70558, 103794,
                                                                       47503, 47653, 72448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 105866, 0, 3,
                                                                       103794, 70684, 103962,
                                                                       47653, 47803, 72658,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 106146, 0, 3,
                                                                       103962, 70810, 104130,
                                                                       47803, 47953, 72868,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 106426, 0, 3,
                                                                       104130, 70936, 104298,
                                                                       47953, 48103, 73078,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 106706, 0, 3,
                                                                       104466, 71608, 104746,
                                                                       48403, 48628, 73918,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 107126, 0, 3,
                                                                       104746, 71818, 105026,
                                                                       48628, 48853, 74233,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 107546, 0, 3,
                                                                       105026, 72028, 105306,
                                                                       48853, 49078, 74548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 107966, 0, 3,
                                                                       105306, 72238, 105586,
                                                                       49078, 49303, 74863,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 108386, 0, 3,
                                                                       105586, 72448, 105866,
                                                                       49303, 49528, 75178,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 108806, 0, 3,
                                                                       105866, 72658, 106146,
                                                                       49528, 49753, 75493,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 109226, 0, 3,
                                                                       106146, 72868, 106426,
                                                                       49753, 49978, 75808,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 109646, 0, 3,
                                                                       106706, 73918, 107126,
                                                                       50428, 50743, 77005,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 110234, 0, 3,
                                                                       107126, 74233, 107546,
                                                                       50743, 51058, 77446,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 110822, 0, 3,
                                                                       107546, 74548, 107966,
                                                                       51058, 51373, 77887,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 111410, 0, 3,
                                                                       107966, 74863, 108386,
                                                                       51373, 51688, 78328,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 111998, 0, 3,
                                                                       108386, 75178, 108806,
                                                                       51688, 52003, 78769,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 112586, 0, 3,
                                                                       108806, 75493, 109226,
                                                                       52003, 52318, 79210,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 113174, 0, 3,
                                                                       109646, 77005, 110234,
                                                                       52948, 53368, 80827,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 113958, 0, 3,
                                                                       110234, 77446, 110822,
                                                                       53368, 53788, 81415,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 114742, 0, 3,
                                                                       110822, 77887, 111410,
                                                                       53788, 54208, 82003,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 115526, 0, 3,
                                                                       111410, 78328, 111998,
                                                                       54208, 54628, 82591,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 116310, 0, 3,
                                                                       111998, 78769, 112586,
                                                                       54628, 55048, 83179,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 117094, 0, 3,
                                                                       113174, 80827, 113958,
                                                                       55888, 56428, 85279,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 118102, 0, 3,
                                                                       113958, 81415, 114742,
                                                                       56428, 56968, 86035,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 119110, 0, 3,
                                                                       114742, 82003, 115526,
                                                                       56968, 57508, 86791,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 120118, 0, 3,
                                                                       115526, 82591, 116310,
                                                                       57508, 58048, 87547,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 121126, 0, 3,
                                                                       117094, 85279, 118102,
                                                                       59128, 59803, 90193,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 122386, 0, 3,
                                                                       118102, 86035, 119110,
                                                                       59803, 60478, 91138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 123646, 0, 3,
                                                                       119110, 86791, 120118,
                                                                       60478, 61153, 92083,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 124906, 0, 3,
                                                                       121126, 90193, 122386,
                                                                       62503, 63328, 95338,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 126446, 0, 3,
                                                                       122386, 91138, 123646,
                                                                       63328, 64153, 96493,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 127986, 0, 3,
                                                                       124906, 95338, 126446,
                                                                       65803, 66793, 100420,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129834, 3, 68773,
                                                                       68794, 101806, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129870, 3, 68794,
                                                                       68815, 101834, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129906, 3, 68815,
                                                                       68836, 101862, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129942, 3, 68836,
                                                                       68857, 101890, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 129978, 3, 68857,
                                                                       68878, 101918, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130014, 3, 68878,
                                                                       68899, 101946, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130050, 3, 68899,
                                                                       68920, 101974, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130086, 3, 68920,
                                                                       68941, 102002, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130122, 3, 68941,
                                                                       68962, 102030, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130158, 3, 68962,
                                                                       68983, 102058, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 130194, 3, 68983,
                                                                       69004, 102086, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130230, 0, 3,
                                                                       129834, 101806, 129870,
                                                                       69046, 69109, 102114,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130338, 0, 3,
                                                                       129870, 101834, 129906,
                                                                       69109, 69172, 102198,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130446, 0, 3,
                                                                       129906, 101862, 129942,
                                                                       69172, 69235, 102282,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130554, 0, 3,
                                                                       129942, 101890, 129978,
                                                                       69235, 69298, 102366,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130662, 0, 3,
                                                                       129978, 101918, 130014,
                                                                       69298, 69361, 102450,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130770, 0, 3,
                                                                       130014, 101946, 130050,
                                                                       69361, 69424, 102534,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130878, 0, 3,
                                                                       130050, 101974, 130086,
                                                                       69424, 69487, 102618,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 130986, 0, 3,
                                                                       130086, 102002, 130122,
                                                                       69487, 69550, 102702,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 131094, 0, 3,
                                                                       130122, 102030, 130158,
                                                                       69550, 69613, 102786,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 131202, 0, 3,
                                                                       130158, 102058, 130194,
                                                                       69613, 69676, 102870,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 131310, 0, 3,
                                                                       130230, 102114, 130338,
                                                                       69802, 69928, 102954,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 131526, 0, 3,
                                                                       130338, 102198, 130446,
                                                                       69928, 70054, 103122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 131742, 0, 3,
                                                                       130446, 102282, 130554,
                                                                       70054, 70180, 103290,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 131958, 0, 3,
                                                                       130554, 102366, 130662,
                                                                       70180, 70306, 103458,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 132174, 0, 3,
                                                                       130662, 102450, 130770,
                                                                       70306, 70432, 103626,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 132390, 0, 3,
                                                                       130770, 102534, 130878,
                                                                       70432, 70558, 103794,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 132606, 0, 3,
                                                                       130878, 102618, 130986,
                                                                       70558, 70684, 103962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 132822, 0, 3,
                                                                       130986, 102702, 131094,
                                                                       70684, 70810, 104130,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 133038, 0, 3,
                                                                       131094, 102786, 131202,
                                                                       70810, 70936, 104298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 133254, 0, 3,
                                                                       131310, 102954, 131526,
                                                                       71188, 71398, 104466,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 133614, 0, 3,
                                                                       131526, 103122, 131742,
                                                                       71398, 71608, 104746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 133974, 0, 3,
                                                                       131742, 103290, 131958,
                                                                       71608, 71818, 105026,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 134334, 0, 3,
                                                                       131958, 103458, 132174,
                                                                       71818, 72028, 105306,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 134694, 0, 3,
                                                                       132174, 103626, 132390,
                                                                       72028, 72238, 105586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 135054, 0, 3,
                                                                       132390, 103794, 132606,
                                                                       72238, 72448, 105866,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 135414, 0, 3,
                                                                       132606, 103962, 132822,
                                                                       72448, 72658, 106146,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 135774, 0, 3,
                                                                       132822, 104130, 133038,
                                                                       72658, 72868, 106426,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 136134, 0, 3,
                                                                       133254, 104466, 133614,
                                                                       73288, 73603, 106706,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 136674, 0, 3,
                                                                       133614, 104746, 133974,
                                                                       73603, 73918, 107126,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 137214, 0, 3,
                                                                       133974, 105026, 134334,
                                                                       73918, 74233, 107546,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 137754, 0, 3,
                                                                       134334, 105306, 134694,
                                                                       74233, 74548, 107966,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 138294, 0, 3,
                                                                       134694, 105586, 135054,
                                                                       74548, 74863, 108386,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 138834, 0, 3,
                                                                       135054, 105866, 135414,
                                                                       74863, 75178, 108806,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 139374, 0, 3,
                                                                       135414, 106146, 135774,
                                                                       75178, 75493, 109226,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 139914, 0, 3,
                                                                       136134, 106706, 136674,
                                                                       76123, 76564, 109646,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 140670, 0, 3,
                                                                       136674, 107126, 137214,
                                                                       76564, 77005, 110234,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 141426, 0, 3,
                                                                       137214, 107546, 137754,
                                                                       77005, 77446, 110822,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 142182, 0, 3,
                                                                       137754, 107966, 138294,
                                                                       77446, 77887, 111410,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 142938, 0, 3,
                                                                       138294, 108386, 138834,
                                                                       77887, 78328, 111998,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 143694, 0, 3,
                                                                       138834, 108806, 139374,
                                                                       78328, 78769, 112586,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 144450, 0, 3,
                                                                       139914, 109646, 140670,
                                                                       79651, 80239, 113174,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 145458, 0, 3,
                                                                       140670, 110234, 141426,
                                                                       80239, 80827, 113958,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 146466, 0, 3,
                                                                       141426, 110822, 142182,
                                                                       80827, 81415, 114742,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 147474, 0, 3,
                                                                       142182, 111410, 142938,
                                                                       81415, 82003, 115526,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 148482, 0, 3,
                                                                       142938, 111998, 143694,
                                                                       82003, 82591, 116310,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 149490, 0, 3,
                                                                       144450, 113174, 145458,
                                                                       83767, 84523, 117094,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 150786, 0, 3,
                                                                       145458, 113958, 146466,
                                                                       84523, 85279, 118102,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 152082, 0, 3,
                                                                       146466, 114742, 147474,
                                                                       85279, 86035, 119110,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 153378, 0, 3,
                                                                       147474, 115526, 148482,
                                                                       86035, 86791, 120118,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 154674, 0, 3,
                                                                       149490, 117094, 150786,
                                                                       88303, 89248, 121126,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 156294, 0, 3,
                                                                       150786, 118102, 152082,
                                                                       89248, 90193, 122386,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 157914, 0, 3,
                                                                       152082, 119110, 153378,
                                                                       90193, 91138, 123646,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 159534, 0, 3,
                                                                       154674, 121126, 156294,
                                                                       93028, 94183, 124906,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 161514, 0, 3,
                                                                       156294, 122386, 157914,
                                                                       94183, 95338, 126446,
                                                                       ncols, gamma, p, q);

                    compute_prim_snk_three_center_electron_repulsion_0(buffer, 163494, 0, 3,
                                                                       159534, 124906, 161514,
                                                                       97648, 99034, 127986,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 165870, 144450, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 167298, 149490, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 169134, 154674, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 171429, 159534, 1980, ncols);

                    simdfunc::contract_primitives(buffer, 174234, 163494, 2376, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 166878, 165870, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 168594, 167298, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 170754, 169134, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 173409, 171429, 55, 1, nmax);

        simdtrf::transform_k_inner(buffer, 176610, 174234, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 177600, 166878, 168594, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 178860, 168594, 170754, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 180480, 170754, 173409, 15, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 182505, 173409, 176610, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 184980, 177600, 178860, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 187500, 178860, 180480, 15, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 190740, 180480, 182505, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 194790, 184980, 187500, 15, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 198990, 187500, 190740, 15, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 204390, 194790, 198990, 15, nmax);

        simdtrf::transform_i_inner(buffer, 210690, 204390, 15, 15, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 210690, 195, nmax);
    }

    for (size_t m = 0; m < 1755; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
