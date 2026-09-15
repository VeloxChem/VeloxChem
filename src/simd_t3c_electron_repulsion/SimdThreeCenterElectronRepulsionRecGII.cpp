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


#include "SimdThreeCenterElectronRepulsionRecGII.hpp"

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
#include "SimdTransferFI.hpp"
#include "SimdTransferFK.hpp"
#include "SimdTransferGI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransferPM.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_gii_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_gii_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 151047, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1521 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 151047, 110404, 8572, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 16,
                                                             ncols, fj, 6, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2725, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2728, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2731, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2734, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2737, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2740, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2743, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2746, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2749, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2752, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2755, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2758, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2761, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2764, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2767, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2770, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2779, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2788, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2797, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2806, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2815, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2824, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2833, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2842, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2851, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2860, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2869, 3, 21, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2878, 3, 22, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2887, 3, 23, 70,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2896, 3, 31, 85,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2914, 3, 34, 91,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2932, 3, 37, 97,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2950, 3, 40, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2968, 3, 43, 109,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2986, 3, 46, 115,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3004, 3, 49, 121,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3022, 3, 52, 127,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3040, 3, 55, 133,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3058, 3, 58, 139,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3076, 3, 61, 145,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3094, 3, 64, 151,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3112, 3, 67, 157,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3130, 3, 85, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3160, 3, 91, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3190, 3, 97, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3220, 3, 103, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3250, 3, 109, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3280, 3, 115, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3310, 3, 121, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3340, 3, 127, 253,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3370, 3, 133, 263,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3400, 3, 139, 273,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3430, 3, 145, 283,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3460, 3, 151, 293,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3490, 3, 183, 333,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3535, 3, 193, 348,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3580, 3, 203, 363,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3625, 3, 213, 378,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3670, 3, 223, 393,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3715, 3, 233, 408,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3760, 3, 243, 423,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3805, 3, 253, 438,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3850, 3, 263, 453,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3895, 3, 273, 468,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3940, 3, 283, 483,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3985, 3, 333, 540,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4048, 3, 348, 561,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4111, 3, 363, 582,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4174, 3, 378, 603,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4237, 3, 393, 624,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4300, 3, 408, 645,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4363, 3, 423, 666,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4426, 3, 438, 687,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4489, 3, 453, 708,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4552, 3, 468, 729,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4615, 3, 540, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4699, 3, 561, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4783, 3, 582, 862,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4867, 3, 603, 890,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4951, 3, 624, 918,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5035, 3, 645, 946,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5119, 3, 666, 974,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5203, 3, 687,
                                                                       1002, ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5287, 3, 708,
                                                                       1030, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5371, 3, 806,
                                                                       1130, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5479, 3, 834,
                                                                       1166, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5587, 3, 862,
                                                                       1202, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5695, 3, 890,
                                                                       1238, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5803, 3, 918,
                                                                       1274, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5911, 3, 946,
                                                                       1310, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6019, 3, 974,
                                                                       1346, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6127, 3, 1002,
                                                                       1382, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6235, 3, 1130,
                                                                       1508, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6370, 3, 1166,
                                                                       1553, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6505, 3, 1202,
                                                                       1598, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6640, 3, 1238,
                                                                       1643, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6775, 3, 1274,
                                                                       1688, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6910, 3, 1310,
                                                                       1733, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7045, 3, 1346,
                                                                       1778, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7180, 3, 1508,
                                                                       1933, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7345, 3, 1553,
                                                                       1988, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7510, 3, 1598,
                                                                       2043, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7675, 3, 1643,
                                                                       2098, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7840, 3, 1688,
                                                                       2153, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8005, 3, 1733,
                                                                       2208, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8170, 3, 1933,
                                                                       2395, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8368, 3, 1988,
                                                                       2461, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8566, 3, 2043,
                                                                       2527, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8764, 3, 2098,
                                                                       2593, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8962, 3, 2153,
                                                                       2659, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9160, 3, 8, 9,
                                                                       2725, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9166, 3, 9, 10,
                                                                       2728, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9172, 3, 10, 11,
                                                                       2731, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9178, 3, 11, 12,
                                                                       2734, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9184, 3, 12, 13,
                                                                       2737, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9190, 3, 13, 14,
                                                                       2740, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9196, 3, 14, 15,
                                                                       2743, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9202, 3, 15, 16,
                                                                       2746, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9208, 3, 16, 17,
                                                                       2749, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9214, 3, 17, 18,
                                                                       2752, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9220, 3, 18, 19,
                                                                       2755, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9226, 3, 19, 20,
                                                                       2758, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9232, 3, 20, 21,
                                                                       2761, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9238, 3, 21, 22,
                                                                       2764, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9244, 3, 22, 23,
                                                                       2767, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9250, 0, 3, 9160,
                                                                       2725, 9166, 2770, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9268, 0, 3, 9166,
                                                                       2728, 9172, 2779, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9286, 0, 3, 9172,
                                                                       2731, 9178, 2788, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9304, 0, 3, 9178,
                                                                       2734, 9184, 2797, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9322, 0, 3, 9184,
                                                                       2737, 9190, 2806, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9340, 0, 3, 9190,
                                                                       2740, 9196, 2815, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9358, 0, 3, 9196,
                                                                       2743, 9202, 2824, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9376, 0, 3, 9202,
                                                                       2746, 9208, 2833, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9394, 0, 3, 9208,
                                                                       2749, 9214, 2842, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9412, 0, 3, 9214,
                                                                       2752, 9220, 2851, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9430, 0, 3, 9220,
                                                                       2755, 9226, 2860, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9448, 0, 3, 9226,
                                                                       2758, 9232, 2869, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9466, 0, 3, 9232,
                                                                       2761, 9238, 2878, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9484, 0, 3, 9238,
                                                                       2764, 9244, 2887, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9502, 0, 3, 9250,
                                                                       2770, 9268, 73, 79, 2896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9538, 0, 3, 9268,
                                                                       2779, 9286, 79, 85, 2914,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9574, 0, 3, 9286,
                                                                       2788, 9304, 85, 91, 2932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9610, 0, 3, 9304,
                                                                       2797, 9322, 91, 97, 2950,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9646, 0, 3, 9322,
                                                                       2806, 9340, 97, 103, 2968,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9682, 0, 3, 9340,
                                                                       2815, 9358, 103, 109,
                                                                       2986, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9718, 0, 3, 9358,
                                                                       2824, 9376, 109, 115,
                                                                       3004, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9754, 0, 3, 9376,
                                                                       2833, 9394, 115, 121,
                                                                       3022, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9790, 0, 3, 9394,
                                                                       2842, 9412, 121, 127,
                                                                       3040, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9826, 0, 3, 9412,
                                                                       2851, 9430, 127, 133,
                                                                       3058, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9862, 0, 3, 9430,
                                                                       2860, 9448, 133, 139,
                                                                       3076, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9898, 0, 3, 9448,
                                                                       2869, 9466, 139, 145,
                                                                       3094, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9934, 0, 3, 9466,
                                                                       2878, 9484, 145, 151,
                                                                       3112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9970, 0, 3, 9502,
                                                                       2896, 9538, 163, 173,
                                                                       3130, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10030, 0, 3, 9538,
                                                                       2914, 9574, 173, 183,
                                                                       3160, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10090, 0, 3, 9574,
                                                                       2932, 9610, 183, 193,
                                                                       3190, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10150, 0, 3, 9610,
                                                                       2950, 9646, 193, 203,
                                                                       3220, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10210, 0, 3, 9646,
                                                                       2968, 9682, 203, 213,
                                                                       3250, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10270, 0, 3, 9682,
                                                                       2986, 9718, 213, 223,
                                                                       3280, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10330, 0, 3, 9718,
                                                                       3004, 9754, 223, 233,
                                                                       3310, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10390, 0, 3, 9754,
                                                                       3022, 9790, 233, 243,
                                                                       3340, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10450, 0, 3, 9790,
                                                                       3040, 9826, 243, 253,
                                                                       3370, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10510, 0, 3, 9826,
                                                                       3058, 9862, 253, 263,
                                                                       3400, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10570, 0, 3, 9862,
                                                                       3076, 9898, 263, 273,
                                                                       3430, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10630, 0, 3, 9898,
                                                                       3094, 9934, 273, 283,
                                                                       3460, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10690, 0, 3, 9970,
                                                                       3130, 10030, 303, 318,
                                                                       3490, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10780, 0, 3,
                                                                       10030, 3160, 10090, 318,
                                                                       333, 3535, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10870, 0, 3,
                                                                       10090, 3190, 10150, 333,
                                                                       348, 3580, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10960, 0, 3,
                                                                       10150, 3220, 10210, 348,
                                                                       363, 3625, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11050, 0, 3,
                                                                       10210, 3250, 10270, 363,
                                                                       378, 3670, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11140, 0, 3,
                                                                       10270, 3280, 10330, 378,
                                                                       393, 3715, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11230, 0, 3,
                                                                       10330, 3310, 10390, 393,
                                                                       408, 3760, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11320, 0, 3,
                                                                       10390, 3340, 10450, 408,
                                                                       423, 3805, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11410, 0, 3,
                                                                       10450, 3370, 10510, 423,
                                                                       438, 3850, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11500, 0, 3,
                                                                       10510, 3400, 10570, 438,
                                                                       453, 3895, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11590, 0, 3,
                                                                       10570, 3430, 10630, 453,
                                                                       468, 3940, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11680, 0, 3,
                                                                       10690, 3490, 10780, 498,
                                                                       519, 3985, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11806, 0, 3,
                                                                       10780, 3535, 10870, 519,
                                                                       540, 4048, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11932, 0, 3,
                                                                       10870, 3580, 10960, 540,
                                                                       561, 4111, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12058, 0, 3,
                                                                       10960, 3625, 11050, 561,
                                                                       582, 4174, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12184, 0, 3,
                                                                       11050, 3670, 11140, 582,
                                                                       603, 4237, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12310, 0, 3,
                                                                       11140, 3715, 11230, 603,
                                                                       624, 4300, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12436, 0, 3,
                                                                       11230, 3760, 11320, 624,
                                                                       645, 4363, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12562, 0, 3,
                                                                       11320, 3805, 11410, 645,
                                                                       666, 4426, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12688, 0, 3,
                                                                       11410, 3850, 11500, 666,
                                                                       687, 4489, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12814, 0, 3,
                                                                       11500, 3895, 11590, 687,
                                                                       708, 4552, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12940, 0, 3,
                                                                       11680, 3985, 11806, 750,
                                                                       778, 4615, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13108, 0, 3,
                                                                       11806, 4048, 11932, 778,
                                                                       806, 4699, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13276, 0, 3,
                                                                       11932, 4111, 12058, 806,
                                                                       834, 4783, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13444, 0, 3,
                                                                       12058, 4174, 12184, 834,
                                                                       862, 4867, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13612, 0, 3,
                                                                       12184, 4237, 12310, 862,
                                                                       890, 4951, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13780, 0, 3,
                                                                       12310, 4300, 12436, 890,
                                                                       918, 5035, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13948, 0, 3,
                                                                       12436, 4363, 12562, 918,
                                                                       946, 5119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14116, 0, 3,
                                                                       12562, 4426, 12688, 946,
                                                                       974, 5203, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14284, 0, 3,
                                                                       12688, 4489, 12814, 974,
                                                                       1002, 5287, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14452, 0, 3,
                                                                       12940, 4615, 13108, 1058,
                                                                       1094, 5371, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14668, 0, 3,
                                                                       13108, 4699, 13276, 1094,
                                                                       1130, 5479, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14884, 0, 3,
                                                                       13276, 4783, 13444, 1130,
                                                                       1166, 5587, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15100, 0, 3,
                                                                       13444, 4867, 13612, 1166,
                                                                       1202, 5695, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15316, 0, 3,
                                                                       13612, 4951, 13780, 1202,
                                                                       1238, 5803, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15532, 0, 3,
                                                                       13780, 5035, 13948, 1238,
                                                                       1274, 5911, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15748, 0, 3,
                                                                       13948, 5119, 14116, 1274,
                                                                       1310, 6019, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15964, 0, 3,
                                                                       14116, 5203, 14284, 1310,
                                                                       1346, 6127, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16180, 0, 3,
                                                                       14452, 5371, 14668, 1418,
                                                                       1463, 6235, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16450, 0, 3,
                                                                       14668, 5479, 14884, 1463,
                                                                       1508, 6370, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16720, 0, 3,
                                                                       14884, 5587, 15100, 1508,
                                                                       1553, 6505, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16990, 0, 3,
                                                                       15100, 5695, 15316, 1553,
                                                                       1598, 6640, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17260, 0, 3,
                                                                       15316, 5803, 15532, 1598,
                                                                       1643, 6775, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17530, 0, 3,
                                                                       15532, 5911, 15748, 1643,
                                                                       1688, 6910, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17800, 0, 3,
                                                                       15748, 6019, 15964, 1688,
                                                                       1733, 7045, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 18070, 0, 3,
                                                                       16180, 6235, 16450, 1823,
                                                                       1878, 7180, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 18400, 0, 3,
                                                                       16450, 6370, 16720, 1878,
                                                                       1933, 7345, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 18730, 0, 3,
                                                                       16720, 6505, 16990, 1933,
                                                                       1988, 7510, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19060, 0, 3,
                                                                       16990, 6640, 17260, 1988,
                                                                       2043, 7675, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19390, 0, 3,
                                                                       17260, 6775, 17530, 2043,
                                                                       2098, 7840, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19720, 0, 3,
                                                                       17530, 6910, 17800, 2098,
                                                                       2153, 8005, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 20050, 0, 3,
                                                                       18070, 7180, 18400, 2263,
                                                                       2329, 8170, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 20446, 0, 3,
                                                                       18400, 7345, 18730, 2329,
                                                                       2395, 8368, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 20842, 0, 3,
                                                                       18730, 7510, 19060, 2395,
                                                                       2461, 8566, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 21238, 0, 3,
                                                                       19060, 7675, 19390, 2461,
                                                                       2527, 8764, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 21634, 0, 3,
                                                                       19390, 7840, 19720, 2527,
                                                                       2593, 8962, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22030, 3, 2725,
                                                                       2728, 9172, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22040, 3, 2728,
                                                                       2731, 9178, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22050, 3, 2731,
                                                                       2734, 9184, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22060, 3, 2734,
                                                                       2737, 9190, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22070, 3, 2737,
                                                                       2740, 9196, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22080, 3, 2740,
                                                                       2743, 9202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22090, 3, 2743,
                                                                       2746, 9208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22100, 3, 2746,
                                                                       2749, 9214, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22110, 3, 2749,
                                                                       2752, 9220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22120, 3, 2752,
                                                                       2755, 9226, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22130, 3, 2755,
                                                                       2758, 9232, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22140, 3, 2758,
                                                                       2761, 9238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 22150, 3, 2761,
                                                                       2764, 9244, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22160, 0, 3,
                                                                       22030, 9172, 22040, 9286,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22190, 0, 3,
                                                                       22040, 9178, 22050, 9304,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22220, 0, 3,
                                                                       22050, 9184, 22060, 9322,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22250, 0, 3,
                                                                       22060, 9190, 22070, 9340,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22280, 0, 3,
                                                                       22070, 9196, 22080, 9358,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22310, 0, 3,
                                                                       22080, 9202, 22090, 9376,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22340, 0, 3,
                                                                       22090, 9208, 22100, 9394,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22370, 0, 3,
                                                                       22100, 9214, 22110, 9412,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22400, 0, 3,
                                                                       22110, 9220, 22120, 9430,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22430, 0, 3,
                                                                       22120, 9226, 22130, 9448,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22460, 0, 3,
                                                                       22130, 9232, 22140, 9466,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22490, 0, 3,
                                                                       22140, 9238, 22150, 9484,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22520, 0, 3,
                                                                       22160, 9286, 22190, 2896,
                                                                       2914, 9574, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22580, 0, 3,
                                                                       22190, 9304, 22220, 2914,
                                                                       2932, 9610, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22640, 0, 3,
                                                                       22220, 9322, 22250, 2932,
                                                                       2950, 9646, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22700, 0, 3,
                                                                       22250, 9340, 22280, 2950,
                                                                       2968, 9682, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22760, 0, 3,
                                                                       22280, 9358, 22310, 2968,
                                                                       2986, 9718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22820, 0, 3,
                                                                       22310, 9376, 22340, 2986,
                                                                       3004, 9754, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22880, 0, 3,
                                                                       22340, 9394, 22370, 3004,
                                                                       3022, 9790, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22940, 0, 3,
                                                                       22370, 9412, 22400, 3022,
                                                                       3040, 9826, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23000, 0, 3,
                                                                       22400, 9430, 22430, 3040,
                                                                       3058, 9862, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23060, 0, 3,
                                                                       22430, 9448, 22460, 3058,
                                                                       3076, 9898, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 23120, 0, 3,
                                                                       22460, 9466, 22490, 3076,
                                                                       3094, 9934, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23180, 0, 3,
                                                                       22520, 9574, 22580, 3130,
                                                                       3160, 10090, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23280, 0, 3,
                                                                       22580, 9610, 22640, 3160,
                                                                       3190, 10150, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23380, 0, 3,
                                                                       22640, 9646, 22700, 3190,
                                                                       3220, 10210, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23480, 0, 3,
                                                                       22700, 9682, 22760, 3220,
                                                                       3250, 10270, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23580, 0, 3,
                                                                       22760, 9718, 22820, 3250,
                                                                       3280, 10330, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23680, 0, 3,
                                                                       22820, 9754, 22880, 3280,
                                                                       3310, 10390, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23780, 0, 3,
                                                                       22880, 9790, 22940, 3310,
                                                                       3340, 10450, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23880, 0, 3,
                                                                       22940, 9826, 23000, 3340,
                                                                       3370, 10510, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23980, 0, 3,
                                                                       23000, 9862, 23060, 3370,
                                                                       3400, 10570, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 24080, 0, 3,
                                                                       23060, 9898, 23120, 3400,
                                                                       3430, 10630, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24180, 0, 3,
                                                                       23180, 10090, 23280, 3490,
                                                                       3535, 10870, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24330, 0, 3,
                                                                       23280, 10150, 23380, 3535,
                                                                       3580, 10960, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24480, 0, 3,
                                                                       23380, 10210, 23480, 3580,
                                                                       3625, 11050, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24630, 0, 3,
                                                                       23480, 10270, 23580, 3625,
                                                                       3670, 11140, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24780, 0, 3,
                                                                       23580, 10330, 23680, 3670,
                                                                       3715, 11230, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24930, 0, 3,
                                                                       23680, 10390, 23780, 3715,
                                                                       3760, 11320, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25080, 0, 3,
                                                                       23780, 10450, 23880, 3760,
                                                                       3805, 11410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25230, 0, 3,
                                                                       23880, 10510, 23980, 3805,
                                                                       3850, 11500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 25380, 0, 3,
                                                                       23980, 10570, 24080, 3850,
                                                                       3895, 11590, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25530, 0, 3,
                                                                       24180, 10870, 24330, 3985,
                                                                       4048, 11932, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25740, 0, 3,
                                                                       24330, 10960, 24480, 4048,
                                                                       4111, 12058, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25950, 0, 3,
                                                                       24480, 11050, 24630, 4111,
                                                                       4174, 12184, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26160, 0, 3,
                                                                       24630, 11140, 24780, 4174,
                                                                       4237, 12310, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26370, 0, 3,
                                                                       24780, 11230, 24930, 4237,
                                                                       4300, 12436, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26580, 0, 3,
                                                                       24930, 11320, 25080, 4300,
                                                                       4363, 12562, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26790, 0, 3,
                                                                       25080, 11410, 25230, 4363,
                                                                       4426, 12688, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27000, 0, 3,
                                                                       25230, 11500, 25380, 4426,
                                                                       4489, 12814, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27210, 0, 3,
                                                                       25530, 11932, 25740, 4615,
                                                                       4699, 13276, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27490, 0, 3,
                                                                       25740, 12058, 25950, 4699,
                                                                       4783, 13444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27770, 0, 3,
                                                                       25950, 12184, 26160, 4783,
                                                                       4867, 13612, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28050, 0, 3,
                                                                       26160, 12310, 26370, 4867,
                                                                       4951, 13780, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28330, 0, 3,
                                                                       26370, 12436, 26580, 4951,
                                                                       5035, 13948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28610, 0, 3,
                                                                       26580, 12562, 26790, 5035,
                                                                       5119, 14116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 28890, 0, 3,
                                                                       26790, 12688, 27000, 5119,
                                                                       5203, 14284, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29170, 0, 3,
                                                                       27210, 13276, 27490, 5371,
                                                                       5479, 14884, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29530, 0, 3,
                                                                       27490, 13444, 27770, 5479,
                                                                       5587, 15100, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29890, 0, 3,
                                                                       27770, 13612, 28050, 5587,
                                                                       5695, 15316, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30250, 0, 3,
                                                                       28050, 13780, 28330, 5695,
                                                                       5803, 15532, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30610, 0, 3,
                                                                       28330, 13948, 28610, 5803,
                                                                       5911, 15748, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 30970, 0, 3,
                                                                       28610, 14116, 28890, 5911,
                                                                       6019, 15964, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 31330, 0, 3,
                                                                       29170, 14884, 29530, 6235,
                                                                       6370, 16720, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 31780, 0, 3,
                                                                       29530, 15100, 29890, 6370,
                                                                       6505, 16990, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 32230, 0, 3,
                                                                       29890, 15316, 30250, 6505,
                                                                       6640, 17260, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 32680, 0, 3,
                                                                       30250, 15532, 30610, 6640,
                                                                       6775, 17530, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 33130, 0, 3,
                                                                       30610, 15748, 30970, 6775,
                                                                       6910, 17800, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 33580, 0, 3,
                                                                       31330, 16720, 31780, 7180,
                                                                       7345, 18730, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 34130, 0, 3,
                                                                       31780, 16990, 32230, 7345,
                                                                       7510, 19060, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 34680, 0, 3,
                                                                       32230, 17260, 32680, 7510,
                                                                       7675, 19390, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 35230, 0, 3,
                                                                       32680, 17530, 33130, 7675,
                                                                       7840, 19720, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 35780, 0, 3,
                                                                       33580, 18730, 34130, 8170,
                                                                       8368, 20842, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 36440, 0, 3,
                                                                       34130, 19060, 34680, 8368,
                                                                       8566, 21238, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 37100, 0, 3,
                                                                       34680, 19390, 35230, 8566,
                                                                       8764, 21634, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37760, 3, 9160,
                                                                       9166, 22030, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37775, 3, 9166,
                                                                       9172, 22040, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37790, 3, 9172,
                                                                       9178, 22050, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37805, 3, 9178,
                                                                       9184, 22060, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37820, 3, 9184,
                                                                       9190, 22070, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37835, 3, 9190,
                                                                       9196, 22080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37850, 3, 9196,
                                                                       9202, 22090, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37865, 3, 9202,
                                                                       9208, 22100, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37880, 3, 9208,
                                                                       9214, 22110, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37895, 3, 9214,
                                                                       9220, 22120, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37910, 3, 9220,
                                                                       9226, 22130, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37925, 3, 9226,
                                                                       9232, 22140, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 37940, 3, 9232,
                                                                       9238, 22150, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 37955, 0, 3,
                                                                       37760, 22030, 37775, 9250,
                                                                       9268, 22160, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38000, 0, 3,
                                                                       37775, 22040, 37790, 9268,
                                                                       9286, 22190, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38045, 0, 3,
                                                                       37790, 22050, 37805, 9286,
                                                                       9304, 22220, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38090, 0, 3,
                                                                       37805, 22060, 37820, 9304,
                                                                       9322, 22250, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38135, 0, 3,
                                                                       37820, 22070, 37835, 9322,
                                                                       9340, 22280, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38180, 0, 3,
                                                                       37835, 22080, 37850, 9340,
                                                                       9358, 22310, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38225, 0, 3,
                                                                       37850, 22090, 37865, 9358,
                                                                       9376, 22340, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38270, 0, 3,
                                                                       37865, 22100, 37880, 9376,
                                                                       9394, 22370, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38315, 0, 3,
                                                                       37880, 22110, 37895, 9394,
                                                                       9412, 22400, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38360, 0, 3,
                                                                       37895, 22120, 37910, 9412,
                                                                       9430, 22430, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38405, 0, 3,
                                                                       37910, 22130, 37925, 9430,
                                                                       9448, 22460, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 38450, 0, 3,
                                                                       37925, 22140, 37940, 9448,
                                                                       9466, 22490, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38495, 0, 3,
                                                                       37955, 22160, 38000, 9502,
                                                                       9538, 22520, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38585, 0, 3,
                                                                       38000, 22190, 38045, 9538,
                                                                       9574, 22580, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38675, 0, 3,
                                                                       38045, 22220, 38090, 9574,
                                                                       9610, 22640, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38765, 0, 3,
                                                                       38090, 22250, 38135, 9610,
                                                                       9646, 22700, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38855, 0, 3,
                                                                       38135, 22280, 38180, 9646,
                                                                       9682, 22760, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 38945, 0, 3,
                                                                       38180, 22310, 38225, 9682,
                                                                       9718, 22820, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39035, 0, 3,
                                                                       38225, 22340, 38270, 9718,
                                                                       9754, 22880, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39125, 0, 3,
                                                                       38270, 22370, 38315, 9754,
                                                                       9790, 22940, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39215, 0, 3,
                                                                       38315, 22400, 38360, 9790,
                                                                       9826, 23000, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39305, 0, 3,
                                                                       38360, 22430, 38405, 9826,
                                                                       9862, 23060, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 39395, 0, 3,
                                                                       38405, 22460, 38450, 9862,
                                                                       9898, 23120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39485, 0, 3,
                                                                       38495, 22520, 38585, 9970,
                                                                       10030, 23180, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39635, 0, 3,
                                                                       38585, 22580, 38675,
                                                                       10030, 10090, 23280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39785, 0, 3,
                                                                       38675, 22640, 38765,
                                                                       10090, 10150, 23380,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 39935, 0, 3,
                                                                       38765, 22700, 38855,
                                                                       10150, 10210, 23480,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40085, 0, 3,
                                                                       38855, 22760, 38945,
                                                                       10210, 10270, 23580,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40235, 0, 3,
                                                                       38945, 22820, 39035,
                                                                       10270, 10330, 23680,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40385, 0, 3,
                                                                       39035, 22880, 39125,
                                                                       10330, 10390, 23780,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40535, 0, 3,
                                                                       39125, 22940, 39215,
                                                                       10390, 10450, 23880,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40685, 0, 3,
                                                                       39215, 23000, 39305,
                                                                       10450, 10510, 23980,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 40835, 0, 3,
                                                                       39305, 23060, 39395,
                                                                       10510, 10570, 24080,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 40985, 0, 3,
                                                                       39485, 23180, 39635,
                                                                       10690, 10780, 24180,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41210, 0, 3,
                                                                       39635, 23280, 39785,
                                                                       10780, 10870, 24330,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41435, 0, 3,
                                                                       39785, 23380, 39935,
                                                                       10870, 10960, 24480,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41660, 0, 3,
                                                                       39935, 23480, 40085,
                                                                       10960, 11050, 24630,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 41885, 0, 3,
                                                                       40085, 23580, 40235,
                                                                       11050, 11140, 24780,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42110, 0, 3,
                                                                       40235, 23680, 40385,
                                                                       11140, 11230, 24930,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42335, 0, 3,
                                                                       40385, 23780, 40535,
                                                                       11230, 11320, 25080,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42560, 0, 3,
                                                                       40535, 23880, 40685,
                                                                       11320, 11410, 25230,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 42785, 0, 3,
                                                                       40685, 23980, 40835,
                                                                       11410, 11500, 25380,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43010, 0, 3,
                                                                       40985, 24180, 41210,
                                                                       11680, 11806, 25530,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43325, 0, 3,
                                                                       41210, 24330, 41435,
                                                                       11806, 11932, 25740,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43640, 0, 3,
                                                                       41435, 24480, 41660,
                                                                       11932, 12058, 25950,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 43955, 0, 3,
                                                                       41660, 24630, 41885,
                                                                       12058, 12184, 26160,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44270, 0, 3,
                                                                       41885, 24780, 42110,
                                                                       12184, 12310, 26370,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44585, 0, 3,
                                                                       42110, 24930, 42335,
                                                                       12310, 12436, 26580,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 44900, 0, 3,
                                                                       42335, 25080, 42560,
                                                                       12436, 12562, 26790,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 45215, 0, 3,
                                                                       42560, 25230, 42785,
                                                                       12562, 12688, 27000,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 45530, 0, 3,
                                                                       43010, 25530, 43325,
                                                                       12940, 13108, 27210,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 45950, 0, 3,
                                                                       43325, 25740, 43640,
                                                                       13108, 13276, 27490,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46370, 0, 3,
                                                                       43640, 25950, 43955,
                                                                       13276, 13444, 27770,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 46790, 0, 3,
                                                                       43955, 26160, 44270,
                                                                       13444, 13612, 28050,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 47210, 0, 3,
                                                                       44270, 26370, 44585,
                                                                       13612, 13780, 28330,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 47630, 0, 3,
                                                                       44585, 26580, 44900,
                                                                       13780, 13948, 28610,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 48050, 0, 3,
                                                                       44900, 26790, 45215,
                                                                       13948, 14116, 28890,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 48470, 0, 3,
                                                                       45530, 27210, 45950,
                                                                       14452, 14668, 29170,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 49010, 0, 3,
                                                                       45950, 27490, 46370,
                                                                       14668, 14884, 29530,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 49550, 0, 3,
                                                                       46370, 27770, 46790,
                                                                       14884, 15100, 29890,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 50090, 0, 3,
                                                                       46790, 28050, 47210,
                                                                       15100, 15316, 30250,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 50630, 0, 3,
                                                                       47210, 28330, 47630,
                                                                       15316, 15532, 30610,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 51170, 0, 3,
                                                                       47630, 28610, 48050,
                                                                       15532, 15748, 30970,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 51710, 0, 3,
                                                                       48470, 29170, 49010,
                                                                       16180, 16450, 31330,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 52385, 0, 3,
                                                                       49010, 29530, 49550,
                                                                       16450, 16720, 31780,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 53060, 0, 3,
                                                                       49550, 29890, 50090,
                                                                       16720, 16990, 32230,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 53735, 0, 3,
                                                                       50090, 30250, 50630,
                                                                       16990, 17260, 32680,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 54410, 0, 3,
                                                                       50630, 30610, 51170,
                                                                       17260, 17530, 33130,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 55085, 0, 3,
                                                                       51710, 31330, 52385,
                                                                       18070, 18400, 33580,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 55910, 0, 3,
                                                                       52385, 31780, 53060,
                                                                       18400, 18730, 34130,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 56735, 0, 3,
                                                                       53060, 32230, 53735,
                                                                       18730, 19060, 34680,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 57560, 0, 3,
                                                                       53735, 32680, 54410,
                                                                       19060, 19390, 35230,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 58385, 0, 3,
                                                                       55085, 33580, 55910,
                                                                       20050, 20446, 35780,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 59375, 0, 3,
                                                                       55910, 34130, 56735,
                                                                       20446, 20842, 36440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 60365, 0, 3,
                                                                       56735, 34680, 57560,
                                                                       20842, 21238, 37100,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61355, 3, 22030,
                                                                       22040, 37790, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61376, 3, 22040,
                                                                       22050, 37805, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61397, 3, 22050,
                                                                       22060, 37820, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61418, 3, 22060,
                                                                       22070, 37835, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61439, 3, 22070,
                                                                       22080, 37850, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61460, 3, 22080,
                                                                       22090, 37865, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61481, 3, 22090,
                                                                       22100, 37880, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61502, 3, 22100,
                                                                       22110, 37895, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61523, 3, 22110,
                                                                       22120, 37910, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61544, 3, 22120,
                                                                       22130, 37925, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 61565, 3, 22130,
                                                                       22140, 37940, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61586, 0, 3,
                                                                       61355, 37790, 61376,
                                                                       22160, 22190, 38045,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61649, 0, 3,
                                                                       61376, 37805, 61397,
                                                                       22190, 22220, 38090,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61712, 0, 3,
                                                                       61397, 37820, 61418,
                                                                       22220, 22250, 38135,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61775, 0, 3,
                                                                       61418, 37835, 61439,
                                                                       22250, 22280, 38180,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61838, 0, 3,
                                                                       61439, 37850, 61460,
                                                                       22280, 22310, 38225,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61901, 0, 3,
                                                                       61460, 37865, 61481,
                                                                       22310, 22340, 38270,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 61964, 0, 3,
                                                                       61481, 37880, 61502,
                                                                       22340, 22370, 38315,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62027, 0, 3,
                                                                       61502, 37895, 61523,
                                                                       22370, 22400, 38360,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62090, 0, 3,
                                                                       61523, 37910, 61544,
                                                                       22400, 22430, 38405,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 62153, 0, 3,
                                                                       61544, 37925, 61565,
                                                                       22430, 22460, 38450,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62216, 0, 3,
                                                                       61586, 38045, 61649,
                                                                       22520, 22580, 38675,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62342, 0, 3,
                                                                       61649, 38090, 61712,
                                                                       22580, 22640, 38765,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62468, 0, 3,
                                                                       61712, 38135, 61775,
                                                                       22640, 22700, 38855,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62594, 0, 3,
                                                                       61775, 38180, 61838,
                                                                       22700, 22760, 38945,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62720, 0, 3,
                                                                       61838, 38225, 61901,
                                                                       22760, 22820, 39035,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62846, 0, 3,
                                                                       61901, 38270, 61964,
                                                                       22820, 22880, 39125,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 62972, 0, 3,
                                                                       61964, 38315, 62027,
                                                                       22880, 22940, 39215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 63098, 0, 3,
                                                                       62027, 38360, 62090,
                                                                       22940, 23000, 39305,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 63224, 0, 3,
                                                                       62090, 38405, 62153,
                                                                       23000, 23060, 39395,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63350, 0, 3,
                                                                       62216, 38675, 62342,
                                                                       23180, 23280, 39785,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63560, 0, 3,
                                                                       62342, 38765, 62468,
                                                                       23280, 23380, 39935,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63770, 0, 3,
                                                                       62468, 38855, 62594,
                                                                       23380, 23480, 40085,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 63980, 0, 3,
                                                                       62594, 38945, 62720,
                                                                       23480, 23580, 40235,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64190, 0, 3,
                                                                       62720, 39035, 62846,
                                                                       23580, 23680, 40385,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64400, 0, 3,
                                                                       62846, 39125, 62972,
                                                                       23680, 23780, 40535,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64610, 0, 3,
                                                                       62972, 39215, 63098,
                                                                       23780, 23880, 40685,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 64820, 0, 3,
                                                                       63098, 39305, 63224,
                                                                       23880, 23980, 40835,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65030, 0, 3,
                                                                       63350, 39785, 63560,
                                                                       24180, 24330, 41435,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65345, 0, 3,
                                                                       63560, 39935, 63770,
                                                                       24330, 24480, 41660,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65660, 0, 3,
                                                                       63770, 40085, 63980,
                                                                       24480, 24630, 41885,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 65975, 0, 3,
                                                                       63980, 40235, 64190,
                                                                       24630, 24780, 42110,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 66290, 0, 3,
                                                                       64190, 40385, 64400,
                                                                       24780, 24930, 42335,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 66605, 0, 3,
                                                                       64400, 40535, 64610,
                                                                       24930, 25080, 42560,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 66920, 0, 3,
                                                                       64610, 40685, 64820,
                                                                       25080, 25230, 42785,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 67235, 0, 3,
                                                                       65030, 41435, 65345,
                                                                       25530, 25740, 43640,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 67676, 0, 3,
                                                                       65345, 41660, 65660,
                                                                       25740, 25950, 43955,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 68117, 0, 3,
                                                                       65660, 41885, 65975,
                                                                       25950, 26160, 44270,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 68558, 0, 3,
                                                                       65975, 42110, 66290,
                                                                       26160, 26370, 44585,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 68999, 0, 3,
                                                                       66290, 42335, 66605,
                                                                       26370, 26580, 44900,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 69440, 0, 3,
                                                                       66605, 42560, 66920,
                                                                       26580, 26790, 45215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 69881, 0, 3,
                                                                       67235, 43640, 67676,
                                                                       27210, 27490, 46370,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 70469, 0, 3,
                                                                       67676, 43955, 68117,
                                                                       27490, 27770, 46790,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 71057, 0, 3,
                                                                       68117, 44270, 68558,
                                                                       27770, 28050, 47210,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 71645, 0, 3,
                                                                       68558, 44585, 68999,
                                                                       28050, 28330, 47630,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 72233, 0, 3,
                                                                       68999, 44900, 69440,
                                                                       28330, 28610, 48050,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 72821, 0, 3,
                                                                       69881, 46370, 70469,
                                                                       29170, 29530, 49550,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 73577, 0, 3,
                                                                       70469, 46790, 71057,
                                                                       29530, 29890, 50090,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 74333, 0, 3,
                                                                       71057, 47210, 71645,
                                                                       29890, 30250, 50630,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 75089, 0, 3,
                                                                       71645, 47630, 72233,
                                                                       30250, 30610, 51170,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 75845, 0, 3,
                                                                       72821, 49550, 73577,
                                                                       31330, 31780, 53060,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 76790, 0, 3,
                                                                       73577, 50090, 74333,
                                                                       31780, 32230, 53735,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 77735, 0, 3,
                                                                       74333, 50630, 75089,
                                                                       32230, 32680, 54410,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 78680, 0, 3,
                                                                       75845, 53060, 76790,
                                                                       33580, 34130, 56735,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 79835, 0, 3,
                                                                       76790, 53735, 77735,
                                                                       34130, 34680, 57560,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 80990, 0, 3,
                                                                       78680, 56735, 79835,
                                                                       35780, 36440, 60365,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82376, 3, 37760,
                                                                       37775, 61355, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82404, 3, 37775,
                                                                       37790, 61376, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82432, 3, 37790,
                                                                       37805, 61397, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82460, 3, 37805,
                                                                       37820, 61418, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82488, 3, 37820,
                                                                       37835, 61439, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82516, 3, 37835,
                                                                       37850, 61460, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82544, 3, 37850,
                                                                       37865, 61481, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82572, 3, 37865,
                                                                       37880, 61502, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82600, 3, 37880,
                                                                       37895, 61523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82628, 3, 37895,
                                                                       37910, 61544, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 82656, 3, 37910,
                                                                       37925, 61565, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 82684, 0, 3,
                                                                       82376, 61355, 82404,
                                                                       37955, 38000, 61586,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 82768, 0, 3,
                                                                       82404, 61376, 82432,
                                                                       38000, 38045, 61649,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 82852, 0, 3,
                                                                       82432, 61397, 82460,
                                                                       38045, 38090, 61712,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 82936, 0, 3,
                                                                       82460, 61418, 82488,
                                                                       38090, 38135, 61775,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83020, 0, 3,
                                                                       82488, 61439, 82516,
                                                                       38135, 38180, 61838,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83104, 0, 3,
                                                                       82516, 61460, 82544,
                                                                       38180, 38225, 61901,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83188, 0, 3,
                                                                       82544, 61481, 82572,
                                                                       38225, 38270, 61964,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83272, 0, 3,
                                                                       82572, 61502, 82600,
                                                                       38270, 38315, 62027,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83356, 0, 3,
                                                                       82600, 61523, 82628,
                                                                       38315, 38360, 62090,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 83440, 0, 3,
                                                                       82628, 61544, 82656,
                                                                       38360, 38405, 62153,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 83524, 0, 3,
                                                                       82684, 61586, 82768,
                                                                       38495, 38585, 62216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 83692, 0, 3,
                                                                       82768, 61649, 82852,
                                                                       38585, 38675, 62342,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 83860, 0, 3,
                                                                       82852, 61712, 82936,
                                                                       38675, 38765, 62468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84028, 0, 3,
                                                                       82936, 61775, 83020,
                                                                       38765, 38855, 62594,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84196, 0, 3,
                                                                       83020, 61838, 83104,
                                                                       38855, 38945, 62720,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84364, 0, 3,
                                                                       83104, 61901, 83188,
                                                                       38945, 39035, 62846,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84532, 0, 3,
                                                                       83188, 61964, 83272,
                                                                       39035, 39125, 62972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84700, 0, 3,
                                                                       83272, 62027, 83356,
                                                                       39125, 39215, 63098,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 84868, 0, 3,
                                                                       83356, 62090, 83440,
                                                                       39215, 39305, 63224,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 85036, 0, 3,
                                                                       83524, 62216, 83692,
                                                                       39485, 39635, 63350,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 85316, 0, 3,
                                                                       83692, 62342, 83860,
                                                                       39635, 39785, 63560,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 85596, 0, 3,
                                                                       83860, 62468, 84028,
                                                                       39785, 39935, 63770,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 85876, 0, 3,
                                                                       84028, 62594, 84196,
                                                                       39935, 40085, 63980,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 86156, 0, 3,
                                                                       84196, 62720, 84364,
                                                                       40085, 40235, 64190,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 86436, 0, 3,
                                                                       84364, 62846, 84532,
                                                                       40235, 40385, 64400,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 86716, 0, 3,
                                                                       84532, 62972, 84700,
                                                                       40385, 40535, 64610,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 86996, 0, 3,
                                                                       84700, 63098, 84868,
                                                                       40535, 40685, 64820,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 87276, 0, 3,
                                                                       85036, 63350, 85316,
                                                                       40985, 41210, 65030,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 87696, 0, 3,
                                                                       85316, 63560, 85596,
                                                                       41210, 41435, 65345,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 88116, 0, 3,
                                                                       85596, 63770, 85876,
                                                                       41435, 41660, 65660,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 88536, 0, 3,
                                                                       85876, 63980, 86156,
                                                                       41660, 41885, 65975,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 88956, 0, 3,
                                                                       86156, 64190, 86436,
                                                                       41885, 42110, 66290,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 89376, 0, 3,
                                                                       86436, 64400, 86716,
                                                                       42110, 42335, 66605,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 89796, 0, 3,
                                                                       86716, 64610, 86996,
                                                                       42335, 42560, 66920,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 90216, 0, 3,
                                                                       87276, 65030, 87696,
                                                                       43010, 43325, 67235,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 90804, 0, 3,
                                                                       87696, 65345, 88116,
                                                                       43325, 43640, 67676,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 91392, 0, 3,
                                                                       88116, 65660, 88536,
                                                                       43640, 43955, 68117,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 91980, 0, 3,
                                                                       88536, 65975, 88956,
                                                                       43955, 44270, 68558,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 92568, 0, 3,
                                                                       88956, 66290, 89376,
                                                                       44270, 44585, 68999,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 93156, 0, 3,
                                                                       89376, 66605, 89796,
                                                                       44585, 44900, 69440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 93744, 0, 3,
                                                                       90216, 67235, 90804,
                                                                       45530, 45950, 69881,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 94528, 0, 3,
                                                                       90804, 67676, 91392,
                                                                       45950, 46370, 70469,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 95312, 0, 3,
                                                                       91392, 68117, 91980,
                                                                       46370, 46790, 71057,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 96096, 0, 3,
                                                                       91980, 68558, 92568,
                                                                       46790, 47210, 71645,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 96880, 0, 3,
                                                                       92568, 68999, 93156,
                                                                       47210, 47630, 72233,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 97664, 0, 3,
                                                                       93744, 69881, 94528,
                                                                       48470, 49010, 72821,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 98672, 0, 3,
                                                                       94528, 70469, 95312,
                                                                       49010, 49550, 73577,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 99680, 0, 3,
                                                                       95312, 71057, 96096,
                                                                       49550, 50090, 74333,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 100688, 0, 3,
                                                                       96096, 71645, 96880,
                                                                       50090, 50630, 75089,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 101696, 0, 3,
                                                                       97664, 72821, 98672,
                                                                       51710, 52385, 75845,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 102956, 0, 3,
                                                                       98672, 73577, 99680,
                                                                       52385, 53060, 76790,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 104216, 0, 3,
                                                                       99680, 74333, 100688,
                                                                       53060, 53735, 77735,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 105476, 0, 3,
                                                                       101696, 75845, 102956,
                                                                       55085, 55910, 78680,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 107016, 0, 3,
                                                                       102956, 76790, 104216,
                                                                       55910, 56735, 79835,
                                                                       ncols, gamma, p, q);

                    compute_prim_sni_three_center_electron_repulsion_0(buffer, 108556, 0, 3,
                                                                       105476, 78680, 107016,
                                                                       58385, 59375, 80990,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 110404, 93744, 784, ncols);

                    simdfunc::contract_primitives(buffer, 111552, 97664, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 113028, 101696, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 114873, 105476, 1540, ncols);

                    simdfunc::contract_primitives(buffer, 117128, 108556, 1848, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 111188, 110404, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 112560, 111552, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 114288, 113028, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 116413, 114873, 55, 1, nmax);

        simdtrf::transform_i_inner(buffer, 118976, 117128, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 119834, 111188, 112560, 13, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 120926, 112560, 114288, 13, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 122330, 114288, 116413, 13, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 124085, 116413, 118976, 13, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 126230, 119834, 120926, 13, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 128414, 120926, 122330, 13, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 131222, 122330, 124085, 13, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 134732, 126230, 128414, 13, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 138372, 128414, 131222, 13, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 143052, 134732, 138372, 13, nmax);

        simdtrf::transform_i_inner(buffer, 148512, 143052, 15, 13, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 148512, 169, nmax);
    }

    for (size_t m = 0; m < 1521; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
