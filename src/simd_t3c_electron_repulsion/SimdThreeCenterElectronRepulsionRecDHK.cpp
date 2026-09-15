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


#include "SimdThreeCenterElectronRepulsionRecDHK.hpp"

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
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_dhk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_dhk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 68210, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 825 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 68210, 58790, 3795, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14},
                                                        ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 61, 0, 3, 8, 9,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 67, 0, 3, 9, 10,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 73, 0, 3, 10, 11,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 79, 0, 3, 11, 12,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 85, 0, 3, 12, 13,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 91, 0, 3, 13, 14,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 97, 0, 3, 14, 15,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 103, 0, 3, 15, 16,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 109, 0, 3, 16, 17,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 115, 0, 3, 17, 18,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 121, 0, 3, 18, 19,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 127, 0, 3, 19, 20,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 22, 25,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 25, 28,
                                                                       67, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 28, 31,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 31, 34,
                                                                       79, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 34, 37,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 37, 40,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 40, 43,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 43, 46,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 46, 49,
                                                                       109, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 49, 52,
                                                                       115, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 233, 0, 3, 52, 55,
                                                                       121, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 61, 67,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 67, 73,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 73, 79,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 79, 85,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 303, 0, 3, 85, 91,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 318, 0, 3, 91, 97,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 333, 0, 3, 97,
                                                                       103, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 348, 0, 3, 103,
                                                                       109, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 363, 0, 3, 109,
                                                                       115, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 378, 0, 3, 115,
                                                                       121, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 393, 0, 3, 133,
                                                                       143, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 414, 0, 3, 143,
                                                                       153, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 435, 0, 3, 153,
                                                                       163, 273, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 456, 0, 3, 163,
                                                                       173, 288, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 477, 0, 3, 173,
                                                                       183, 303, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 498, 0, 3, 183,
                                                                       193, 318, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 519, 0, 3, 193,
                                                                       203, 333, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 540, 0, 3, 203,
                                                                       213, 348, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 561, 0, 3, 213,
                                                                       223, 363, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 582, 0, 3, 243,
                                                                       258, 393, 414, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 610, 0, 3, 258,
                                                                       273, 414, 435, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 638, 0, 3, 273,
                                                                       288, 435, 456, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 666, 0, 3, 288,
                                                                       303, 456, 477, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 694, 0, 3, 303,
                                                                       318, 477, 498, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 722, 0, 3, 318,
                                                                       333, 498, 519, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 750, 0, 3, 333,
                                                                       348, 519, 540, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 778, 0, 3, 348,
                                                                       363, 540, 561, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 806, 0, 3, 393,
                                                                       414, 582, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 842, 0, 3, 414,
                                                                       435, 610, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 878, 0, 3, 435,
                                                                       456, 638, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 914, 0, 3, 456,
                                                                       477, 666, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 950, 0, 3, 477,
                                                                       498, 694, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 986, 0, 3, 498,
                                                                       519, 722, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1022, 0, 3, 519,
                                                                       540, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1058, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1061, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1064, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1067, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1070, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1073, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1076, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1079, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1082, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1085, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1088, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1091, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1094, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1097, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1100, 3, 10, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1109, 3, 11, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1118, 3, 12, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1127, 3, 13, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1136, 3, 14, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1145, 3, 15, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1154, 3, 16, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1163, 3, 17, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1172, 3, 18, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1181, 3, 19, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1190, 3, 20, 58,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1199, 3, 22, 61,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1217, 3, 25, 67,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1235, 3, 28, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1253, 3, 31, 79,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1271, 3, 34, 85,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1289, 3, 37, 91,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1307, 3, 40, 97,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1325, 3, 43, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1343, 3, 46, 109,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1361, 3, 49, 115,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1379, 3, 52, 121,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1397, 3, 55, 127,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1415, 3, 61, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1445, 3, 67, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1475, 3, 73, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1505, 3, 79, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1535, 3, 85, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1565, 3, 91, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1595, 3, 97, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1625, 3, 103, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1655, 3, 109, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1685, 3, 115, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1715, 3, 121, 233,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1745, 3, 133, 243,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1790, 3, 143, 258,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1835, 3, 153, 273,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1880, 3, 163, 288,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1925, 3, 173, 303,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1970, 3, 183, 318,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2015, 3, 193, 333,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2060, 3, 203, 348,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2105, 3, 213, 363,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2150, 3, 223, 378,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2195, 3, 243, 393,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2258, 3, 258, 414,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2321, 3, 273, 435,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2384, 3, 288, 456,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2447, 3, 303, 477,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2510, 3, 318, 498,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2573, 3, 333, 519,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2636, 3, 348, 540,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2699, 3, 363, 561,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2762, 3, 393, 582,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2846, 3, 414, 610,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 2930, 3, 435, 638,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3014, 3, 456, 666,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3098, 3, 477, 694,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3182, 3, 498, 722,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3266, 3, 519, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3350, 3, 540, 778,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3434, 3, 582, 806,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3542, 3, 610, 842,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3650, 3, 638, 878,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3758, 3, 666, 914,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3866, 3, 694, 950,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3974, 3, 722, 986,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4082, 3, 750,
                                                                       1022, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4190, 3, 8, 9,
                                                                       1064, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4196, 3, 9, 10,
                                                                       1067, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4202, 3, 10, 11,
                                                                       1070, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4208, 3, 11, 12,
                                                                       1073, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4214, 3, 12, 13,
                                                                       1076, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4220, 3, 13, 14,
                                                                       1079, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4226, 3, 14, 15,
                                                                       1082, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4232, 3, 15, 16,
                                                                       1085, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4238, 3, 16, 17,
                                                                       1088, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4244, 3, 17, 18,
                                                                       1091, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4250, 3, 18, 19,
                                                                       1094, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4256, 3, 19, 20,
                                                                       1097, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4262, 0, 3, 4190,
                                                                       1064, 4196, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4280, 0, 3, 4196,
                                                                       1067, 4202, 1109, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4298, 0, 3, 4202,
                                                                       1070, 4208, 1118, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4316, 0, 3, 4208,
                                                                       1073, 4214, 1127, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4334, 0, 3, 4214,
                                                                       1076, 4220, 1136, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4352, 0, 3, 4220,
                                                                       1079, 4226, 1145, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4370, 0, 3, 4226,
                                                                       1082, 4232, 1154, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4388, 0, 3, 4232,
                                                                       1085, 4238, 1163, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4406, 0, 3, 4238,
                                                                       1088, 4244, 1172, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4424, 0, 3, 4244,
                                                                       1091, 4250, 1181, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4442, 0, 3, 4250,
                                                                       1094, 4256, 1190, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4460, 0, 3, 4262,
                                                                       1100, 4280, 61, 67, 1235,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4496, 0, 3, 4280,
                                                                       1109, 4298, 67, 73, 1253,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4532, 0, 3, 4298,
                                                                       1118, 4316, 73, 79, 1271,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4568, 0, 3, 4316,
                                                                       1127, 4334, 79, 85, 1289,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4604, 0, 3, 4334,
                                                                       1136, 4352, 85, 91, 1307,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4640, 0, 3, 4352,
                                                                       1145, 4370, 91, 97, 1325,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4676, 0, 3, 4370,
                                                                       1154, 4388, 97, 103, 1343,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4712, 0, 3, 4388,
                                                                       1163, 4406, 103, 109,
                                                                       1361, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4748, 0, 3, 4406,
                                                                       1172, 4424, 109, 115,
                                                                       1379, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4784, 0, 3, 4424,
                                                                       1181, 4442, 115, 121,
                                                                       1397, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4820, 0, 3, 4460,
                                                                       1235, 4496, 133, 143,
                                                                       1475, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4880, 0, 3, 4496,
                                                                       1253, 4532, 143, 153,
                                                                       1505, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4940, 0, 3, 4532,
                                                                       1271, 4568, 153, 163,
                                                                       1535, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5000, 0, 3, 4568,
                                                                       1289, 4604, 163, 173,
                                                                       1565, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5060, 0, 3, 4604,
                                                                       1307, 4640, 173, 183,
                                                                       1595, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5120, 0, 3, 4640,
                                                                       1325, 4676, 183, 193,
                                                                       1625, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5180, 0, 3, 4676,
                                                                       1343, 4712, 193, 203,
                                                                       1655, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5240, 0, 3, 4712,
                                                                       1361, 4748, 203, 213,
                                                                       1685, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5300, 0, 3, 4748,
                                                                       1379, 4784, 213, 223,
                                                                       1715, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5360, 0, 3, 4820,
                                                                       1475, 4880, 243, 258,
                                                                       1835, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5450, 0, 3, 4880,
                                                                       1505, 4940, 258, 273,
                                                                       1880, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5540, 0, 3, 4940,
                                                                       1535, 5000, 273, 288,
                                                                       1925, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5630, 0, 3, 5000,
                                                                       1565, 5060, 288, 303,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5720, 0, 3, 5060,
                                                                       1595, 5120, 303, 318,
                                                                       2015, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5810, 0, 3, 5120,
                                                                       1625, 5180, 318, 333,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5900, 0, 3, 5180,
                                                                       1655, 5240, 333, 348,
                                                                       2105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5990, 0, 3, 5240,
                                                                       1685, 5300, 348, 363,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6080, 0, 3, 5360,
                                                                       1835, 5450, 393, 414,
                                                                       2321, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6206, 0, 3, 5450,
                                                                       1880, 5540, 414, 435,
                                                                       2384, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6332, 0, 3, 5540,
                                                                       1925, 5630, 435, 456,
                                                                       2447, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6458, 0, 3, 5630,
                                                                       1970, 5720, 456, 477,
                                                                       2510, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6584, 0, 3, 5720,
                                                                       2015, 5810, 477, 498,
                                                                       2573, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6710, 0, 3, 5810,
                                                                       2060, 5900, 498, 519,
                                                                       2636, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6836, 0, 3, 5900,
                                                                       2105, 5990, 519, 540,
                                                                       2699, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 6962, 0, 3, 6080,
                                                                       2321, 6206, 582, 610,
                                                                       2930, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7130, 0, 3, 6206,
                                                                       2384, 6332, 610, 638,
                                                                       3014, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7298, 0, 3, 6332,
                                                                       2447, 6458, 638, 666,
                                                                       3098, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7466, 0, 3, 6458,
                                                                       2510, 6584, 666, 694,
                                                                       3182, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7634, 0, 3, 6584,
                                                                       2573, 6710, 694, 722,
                                                                       3266, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7802, 0, 3, 6710,
                                                                       2636, 6836, 722, 750,
                                                                       3350, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 7970, 0, 3, 6962,
                                                                       2930, 7130, 806, 842,
                                                                       3650, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 8186, 0, 3, 7130,
                                                                       3014, 7298, 842, 878,
                                                                       3758, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 8402, 0, 3, 7298,
                                                                       3098, 7466, 878, 914,
                                                                       3866, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 8618, 0, 3, 7466,
                                                                       3182, 7634, 914, 950,
                                                                       3974, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 8834, 0, 3, 7634,
                                                                       3266, 7802, 950, 986,
                                                                       4082, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9050, 3, 1058,
                                                                       1061, 4190, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9060, 3, 1061,
                                                                       1064, 4196, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9070, 3, 1064,
                                                                       1067, 4202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9080, 3, 1067,
                                                                       1070, 4208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9090, 3, 1070,
                                                                       1073, 4214, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9100, 3, 1073,
                                                                       1076, 4220, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9110, 3, 1076,
                                                                       1079, 4226, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9120, 3, 1079,
                                                                       1082, 4232, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9130, 3, 1082,
                                                                       1085, 4238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9140, 3, 1085,
                                                                       1088, 4244, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9150, 3, 1088,
                                                                       1091, 4250, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9160, 3, 1091,
                                                                       1094, 4256, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9170, 0, 3, 9050,
                                                                       4190, 9060, 4262, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9200, 0, 3, 9060,
                                                                       4196, 9070, 4280, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9230, 0, 3, 9070,
                                                                       4202, 9080, 4298, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9260, 0, 3, 9080,
                                                                       4208, 9090, 4316, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9290, 0, 3, 9090,
                                                                       4214, 9100, 4334, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9320, 0, 3, 9100,
                                                                       4220, 9110, 4352, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9350, 0, 3, 9110,
                                                                       4226, 9120, 4370, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9380, 0, 3, 9120,
                                                                       4232, 9130, 4388, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9410, 0, 3, 9130,
                                                                       4238, 9140, 4406, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9440, 0, 3, 9140,
                                                                       4244, 9150, 4424, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9470, 0, 3, 9150,
                                                                       4250, 9160, 4442, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9500, 0, 3, 9170,
                                                                       4262, 9200, 1199, 1217,
                                                                       4460, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9560, 0, 3, 9200,
                                                                       4280, 9230, 1217, 1235,
                                                                       4496, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9620, 0, 3, 9230,
                                                                       4298, 9260, 1235, 1253,
                                                                       4532, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9680, 0, 3, 9260,
                                                                       4316, 9290, 1253, 1271,
                                                                       4568, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9740, 0, 3, 9290,
                                                                       4334, 9320, 1271, 1289,
                                                                       4604, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9800, 0, 3, 9320,
                                                                       4352, 9350, 1289, 1307,
                                                                       4640, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9860, 0, 3, 9350,
                                                                       4370, 9380, 1307, 1325,
                                                                       4676, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9920, 0, 3, 9380,
                                                                       4388, 9410, 1325, 1343,
                                                                       4712, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9980, 0, 3, 9410,
                                                                       4406, 9440, 1343, 1361,
                                                                       4748, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10040, 0, 3, 9440,
                                                                       4424, 9470, 1361, 1379,
                                                                       4784, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10100, 0, 3, 9500,
                                                                       4460, 9560, 1415, 1445,
                                                                       4820, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10200, 0, 3, 9560,
                                                                       4496, 9620, 1445, 1475,
                                                                       4880, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10300, 0, 3, 9620,
                                                                       4532, 9680, 1475, 1505,
                                                                       4940, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10400, 0, 3, 9680,
                                                                       4568, 9740, 1505, 1535,
                                                                       5000, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10500, 0, 3, 9740,
                                                                       4604, 9800, 1535, 1565,
                                                                       5060, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10600, 0, 3, 9800,
                                                                       4640, 9860, 1565, 1595,
                                                                       5120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10700, 0, 3, 9860,
                                                                       4676, 9920, 1595, 1625,
                                                                       5180, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10800, 0, 3, 9920,
                                                                       4712, 9980, 1625, 1655,
                                                                       5240, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 10900, 0, 3, 9980,
                                                                       4748, 10040, 1655, 1685,
                                                                       5300, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11000, 0, 3,
                                                                       10100, 4820, 10200, 1745,
                                                                       1790, 5360, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11150, 0, 3,
                                                                       10200, 4880, 10300, 1790,
                                                                       1835, 5450, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11300, 0, 3,
                                                                       10300, 4940, 10400, 1835,
                                                                       1880, 5540, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11450, 0, 3,
                                                                       10400, 5000, 10500, 1880,
                                                                       1925, 5630, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11600, 0, 3,
                                                                       10500, 5060, 10600, 1925,
                                                                       1970, 5720, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11750, 0, 3,
                                                                       10600, 5120, 10700, 1970,
                                                                       2015, 5810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 11900, 0, 3,
                                                                       10700, 5180, 10800, 2015,
                                                                       2060, 5900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12050, 0, 3,
                                                                       10800, 5240, 10900, 2060,
                                                                       2105, 5990, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12200, 0, 3,
                                                                       11000, 5360, 11150, 2195,
                                                                       2258, 6080, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12410, 0, 3,
                                                                       11150, 5450, 11300, 2258,
                                                                       2321, 6206, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12620, 0, 3,
                                                                       11300, 5540, 11450, 2321,
                                                                       2384, 6332, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 12830, 0, 3,
                                                                       11450, 5630, 11600, 2384,
                                                                       2447, 6458, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13040, 0, 3,
                                                                       11600, 5720, 11750, 2447,
                                                                       2510, 6584, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13250, 0, 3,
                                                                       11750, 5810, 11900, 2510,
                                                                       2573, 6710, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13460, 0, 3,
                                                                       11900, 5900, 12050, 2573,
                                                                       2636, 6836, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 13670, 0, 3,
                                                                       12200, 6080, 12410, 2762,
                                                                       2846, 6962, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 13950, 0, 3,
                                                                       12410, 6206, 12620, 2846,
                                                                       2930, 7130, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 14230, 0, 3,
                                                                       12620, 6332, 12830, 2930,
                                                                       3014, 7298, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 14510, 0, 3,
                                                                       12830, 6458, 13040, 3014,
                                                                       3098, 7466, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 14790, 0, 3,
                                                                       13040, 6584, 13250, 3098,
                                                                       3182, 7634, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 15070, 0, 3,
                                                                       13250, 6710, 13460, 3182,
                                                                       3266, 7802, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 15350, 0, 3,
                                                                       13670, 6962, 13950, 3434,
                                                                       3542, 7970, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 15710, 0, 3,
                                                                       13950, 7130, 14230, 3542,
                                                                       3650, 8186, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 16070, 0, 3,
                                                                       14230, 7298, 14510, 3650,
                                                                       3758, 8402, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 16430, 0, 3,
                                                                       14510, 7466, 14790, 3758,
                                                                       3866, 8618, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 16790, 0, 3,
                                                                       14790, 7634, 15070, 3866,
                                                                       3974, 8834, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17150, 3, 4190,
                                                                       4196, 9070, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17165, 3, 4196,
                                                                       4202, 9080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17180, 3, 4202,
                                                                       4208, 9090, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17195, 3, 4208,
                                                                       4214, 9100, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17210, 3, 4214,
                                                                       4220, 9110, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17225, 3, 4220,
                                                                       4226, 9120, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17240, 3, 4226,
                                                                       4232, 9130, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17255, 3, 4232,
                                                                       4238, 9140, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17270, 3, 4238,
                                                                       4244, 9150, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 17285, 3, 4244,
                                                                       4250, 9160, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17300, 0, 3,
                                                                       17150, 9070, 17165, 4262,
                                                                       4280, 9230, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17345, 0, 3,
                                                                       17165, 9080, 17180, 4280,
                                                                       4298, 9260, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17390, 0, 3,
                                                                       17180, 9090, 17195, 4298,
                                                                       4316, 9290, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17435, 0, 3,
                                                                       17195, 9100, 17210, 4316,
                                                                       4334, 9320, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17480, 0, 3,
                                                                       17210, 9110, 17225, 4334,
                                                                       4352, 9350, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17525, 0, 3,
                                                                       17225, 9120, 17240, 4352,
                                                                       4370, 9380, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17570, 0, 3,
                                                                       17240, 9130, 17255, 4370,
                                                                       4388, 9410, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17615, 0, 3,
                                                                       17255, 9140, 17270, 4388,
                                                                       4406, 9440, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17660, 0, 3,
                                                                       17270, 9150, 17285, 4406,
                                                                       4424, 9470, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17705, 0, 3,
                                                                       17300, 9230, 17345, 4460,
                                                                       4496, 9620, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17795, 0, 3,
                                                                       17345, 9260, 17390, 4496,
                                                                       4532, 9680, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17885, 0, 3,
                                                                       17390, 9290, 17435, 4532,
                                                                       4568, 9740, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17975, 0, 3,
                                                                       17435, 9320, 17480, 4568,
                                                                       4604, 9800, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18065, 0, 3,
                                                                       17480, 9350, 17525, 4604,
                                                                       4640, 9860, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18155, 0, 3,
                                                                       17525, 9380, 17570, 4640,
                                                                       4676, 9920, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18245, 0, 3,
                                                                       17570, 9410, 17615, 4676,
                                                                       4712, 9980, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 18335, 0, 3,
                                                                       17615, 9440, 17660, 4712,
                                                                       4748, 10040, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18425, 0, 3,
                                                                       17705, 9620, 17795, 4820,
                                                                       4880, 10300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18575, 0, 3,
                                                                       17795, 9680, 17885, 4880,
                                                                       4940, 10400, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18725, 0, 3,
                                                                       17885, 9740, 17975, 4940,
                                                                       5000, 10500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18875, 0, 3,
                                                                       17975, 9800, 18065, 5000,
                                                                       5060, 10600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19025, 0, 3,
                                                                       18065, 9860, 18155, 5060,
                                                                       5120, 10700, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19175, 0, 3,
                                                                       18155, 9920, 18245, 5120,
                                                                       5180, 10800, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19325, 0, 3,
                                                                       18245, 9980, 18335, 5180,
                                                                       5240, 10900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 19475, 0, 3,
                                                                       18425, 10300, 18575, 5360,
                                                                       5450, 11300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 19700, 0, 3,
                                                                       18575, 10400, 18725, 5450,
                                                                       5540, 11450, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 19925, 0, 3,
                                                                       18725, 10500, 18875, 5540,
                                                                       5630, 11600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 20150, 0, 3,
                                                                       18875, 10600, 19025, 5630,
                                                                       5720, 11750, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 20375, 0, 3,
                                                                       19025, 10700, 19175, 5720,
                                                                       5810, 11900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 20600, 0, 3,
                                                                       19175, 10800, 19325, 5810,
                                                                       5900, 12050, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 20825, 0, 3,
                                                                       19475, 11300, 19700, 6080,
                                                                       6206, 12620, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 21140, 0, 3,
                                                                       19700, 11450, 19925, 6206,
                                                                       6332, 12830, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 21455, 0, 3,
                                                                       19925, 11600, 20150, 6332,
                                                                       6458, 13040, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 21770, 0, 3,
                                                                       20150, 11750, 20375, 6458,
                                                                       6584, 13250, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 22085, 0, 3,
                                                                       20375, 11900, 20600, 6584,
                                                                       6710, 13460, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 22400, 0, 3,
                                                                       20825, 12620, 21140, 6962,
                                                                       7130, 14230, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 22820, 0, 3,
                                                                       21140, 12830, 21455, 7130,
                                                                       7298, 14510, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 23240, 0, 3,
                                                                       21455, 13040, 21770, 7298,
                                                                       7466, 14790, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 23660, 0, 3,
                                                                       21770, 13250, 22085, 7466,
                                                                       7634, 15070, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 24080, 0, 3,
                                                                       22400, 14230, 22820, 7970,
                                                                       8186, 16070, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 24620, 0, 3,
                                                                       22820, 14510, 23240, 8186,
                                                                       8402, 16430, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 25160, 0, 3,
                                                                       23240, 14790, 23660, 8402,
                                                                       8618, 16790, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25700, 3, 9050,
                                                                       9060, 17150, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25721, 3, 9060,
                                                                       9070, 17165, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25742, 3, 9070,
                                                                       9080, 17180, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25763, 3, 9080,
                                                                       9090, 17195, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25784, 3, 9090,
                                                                       9100, 17210, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25805, 3, 9100,
                                                                       9110, 17225, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25826, 3, 9110,
                                                                       9120, 17240, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25847, 3, 9120,
                                                                       9130, 17255, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25868, 3, 9130,
                                                                       9140, 17270, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25889, 3, 9140,
                                                                       9150, 17285, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 25910, 0, 3,
                                                                       25700, 17150, 25721, 9170,
                                                                       9200, 17300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 25973, 0, 3,
                                                                       25721, 17165, 25742, 9200,
                                                                       9230, 17345, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26036, 0, 3,
                                                                       25742, 17180, 25763, 9230,
                                                                       9260, 17390, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26099, 0, 3,
                                                                       25763, 17195, 25784, 9260,
                                                                       9290, 17435, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26162, 0, 3,
                                                                       25784, 17210, 25805, 9290,
                                                                       9320, 17480, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26225, 0, 3,
                                                                       25805, 17225, 25826, 9320,
                                                                       9350, 17525, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26288, 0, 3,
                                                                       25826, 17240, 25847, 9350,
                                                                       9380, 17570, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26351, 0, 3,
                                                                       25847, 17255, 25868, 9380,
                                                                       9410, 17615, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26414, 0, 3,
                                                                       25868, 17270, 25889, 9410,
                                                                       9440, 17660, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26477, 0, 3,
                                                                       25910, 17300, 25973, 9500,
                                                                       9560, 17705, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26603, 0, 3,
                                                                       25973, 17345, 26036, 9560,
                                                                       9620, 17795, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26729, 0, 3,
                                                                       26036, 17390, 26099, 9620,
                                                                       9680, 17885, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26855, 0, 3,
                                                                       26099, 17435, 26162, 9680,
                                                                       9740, 17975, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26981, 0, 3,
                                                                       26162, 17480, 26225, 9740,
                                                                       9800, 18065, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 27107, 0, 3,
                                                                       26225, 17525, 26288, 9800,
                                                                       9860, 18155, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 27233, 0, 3,
                                                                       26288, 17570, 26351, 9860,
                                                                       9920, 18245, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 27359, 0, 3,
                                                                       26351, 17615, 26414, 9920,
                                                                       9980, 18335, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27485, 0, 3,
                                                                       26477, 17705, 26603,
                                                                       10100, 10200, 18425,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27695, 0, 3,
                                                                       26603, 17795, 26729,
                                                                       10200, 10300, 18575,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27905, 0, 3,
                                                                       26729, 17885, 26855,
                                                                       10300, 10400, 18725,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 28115, 0, 3,
                                                                       26855, 17975, 26981,
                                                                       10400, 10500, 18875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 28325, 0, 3,
                                                                       26981, 18065, 27107,
                                                                       10500, 10600, 19025,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 28535, 0, 3,
                                                                       27107, 18155, 27233,
                                                                       10600, 10700, 19175,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 28745, 0, 3,
                                                                       27233, 18245, 27359,
                                                                       10700, 10800, 19325,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 28955, 0, 3,
                                                                       27485, 18425, 27695,
                                                                       11000, 11150, 19475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 29270, 0, 3,
                                                                       27695, 18575, 27905,
                                                                       11150, 11300, 19700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 29585, 0, 3,
                                                                       27905, 18725, 28115,
                                                                       11300, 11450, 19925,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 29900, 0, 3,
                                                                       28115, 18875, 28325,
                                                                       11450, 11600, 20150,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 30215, 0, 3,
                                                                       28325, 19025, 28535,
                                                                       11600, 11750, 20375,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 30530, 0, 3,
                                                                       28535, 19175, 28745,
                                                                       11750, 11900, 20600,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 30845, 0, 3,
                                                                       28955, 19475, 29270,
                                                                       12200, 12410, 20825,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 31286, 0, 3,
                                                                       29270, 19700, 29585,
                                                                       12410, 12620, 21140,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 31727, 0, 3,
                                                                       29585, 19925, 29900,
                                                                       12620, 12830, 21455,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 32168, 0, 3,
                                                                       29900, 20150, 30215,
                                                                       12830, 13040, 21770,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 32609, 0, 3,
                                                                       30215, 20375, 30530,
                                                                       13040, 13250, 22085,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 33050, 0, 3,
                                                                       30845, 20825, 31286,
                                                                       13670, 13950, 22400,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 33638, 0, 3,
                                                                       31286, 21140, 31727,
                                                                       13950, 14230, 22820,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 34226, 0, 3,
                                                                       31727, 21455, 32168,
                                                                       14230, 14510, 23240,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 34814, 0, 3,
                                                                       32168, 21770, 32609,
                                                                       14510, 14790, 23660,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 35402, 0, 3,
                                                                       33050, 22400, 33638,
                                                                       15350, 15710, 24080,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 36158, 0, 3,
                                                                       33638, 22820, 34226,
                                                                       15710, 16070, 24620,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 36914, 0, 3,
                                                                       34226, 23240, 34814,
                                                                       16070, 16430, 25160,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37670, 3, 17150,
                                                                       17165, 25742, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37698, 3, 17165,
                                                                       17180, 25763, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37726, 3, 17180,
                                                                       17195, 25784, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37754, 3, 17195,
                                                                       17210, 25805, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37782, 3, 17210,
                                                                       17225, 25826, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37810, 3, 17225,
                                                                       17240, 25847, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37838, 3, 17240,
                                                                       17255, 25868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 37866, 3, 17255,
                                                                       17270, 25889, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 37894, 0, 3,
                                                                       37670, 25742, 37698,
                                                                       17300, 17345, 26036,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 37978, 0, 3,
                                                                       37698, 25763, 37726,
                                                                       17345, 17390, 26099,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38062, 0, 3,
                                                                       37726, 25784, 37754,
                                                                       17390, 17435, 26162,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38146, 0, 3,
                                                                       37754, 25805, 37782,
                                                                       17435, 17480, 26225,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38230, 0, 3,
                                                                       37782, 25826, 37810,
                                                                       17480, 17525, 26288,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38314, 0, 3,
                                                                       37810, 25847, 37838,
                                                                       17525, 17570, 26351,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 38398, 0, 3,
                                                                       37838, 25868, 37866,
                                                                       17570, 17615, 26414,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 38482, 0, 3,
                                                                       37894, 26036, 37978,
                                                                       17705, 17795, 26729,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 38650, 0, 3,
                                                                       37978, 26099, 38062,
                                                                       17795, 17885, 26855,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 38818, 0, 3,
                                                                       38062, 26162, 38146,
                                                                       17885, 17975, 26981,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 38986, 0, 3,
                                                                       38146, 26225, 38230,
                                                                       17975, 18065, 27107,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 39154, 0, 3,
                                                                       38230, 26288, 38314,
                                                                       18065, 18155, 27233,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 39322, 0, 3,
                                                                       38314, 26351, 38398,
                                                                       18155, 18245, 27359,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 39490, 0, 3,
                                                                       38482, 26729, 38650,
                                                                       18425, 18575, 27905,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 39770, 0, 3,
                                                                       38650, 26855, 38818,
                                                                       18575, 18725, 28115,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 40050, 0, 3,
                                                                       38818, 26981, 38986,
                                                                       18725, 18875, 28325,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 40330, 0, 3,
                                                                       38986, 27107, 39154,
                                                                       18875, 19025, 28535,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 40610, 0, 3,
                                                                       39154, 27233, 39322,
                                                                       19025, 19175, 28745,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 40890, 0, 3,
                                                                       39490, 27905, 39770,
                                                                       19475, 19700, 29585,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 41310, 0, 3,
                                                                       39770, 28115, 40050,
                                                                       19700, 19925, 29900,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 41730, 0, 3,
                                                                       40050, 28325, 40330,
                                                                       19925, 20150, 30215,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 42150, 0, 3,
                                                                       40330, 28535, 40610,
                                                                       20150, 20375, 30530,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 42570, 0, 3,
                                                                       40890, 29585, 41310,
                                                                       20825, 21140, 31727,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 43158, 0, 3,
                                                                       41310, 29900, 41730,
                                                                       21140, 21455, 32168,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 43746, 0, 3,
                                                                       41730, 30215, 42150,
                                                                       21455, 21770, 32609,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 44334, 0, 3,
                                                                       42570, 31727, 43158,
                                                                       22400, 22820, 34226,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 45118, 0, 3,
                                                                       43158, 32168, 43746,
                                                                       22820, 23240, 34814,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 45902, 0, 3,
                                                                       44334, 34226, 45118,
                                                                       24080, 24620, 36914,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 46910, 3, 25700,
                                                                       25721, 37670, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 46946, 3, 25721,
                                                                       25742, 37698, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 46982, 3, 25742,
                                                                       25763, 37726, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47018, 3, 25763,
                                                                       25784, 37754, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47054, 3, 25784,
                                                                       25805, 37782, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47090, 3, 25805,
                                                                       25826, 37810, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47126, 3, 25826,
                                                                       25847, 37838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 47162, 3, 25847,
                                                                       25868, 37866, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 47198, 0, 3,
                                                                       46910, 37670, 46946,
                                                                       25910, 25973, 37894,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 47306, 0, 3,
                                                                       46946, 37698, 46982,
                                                                       25973, 26036, 37978,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 47414, 0, 3,
                                                                       46982, 37726, 47018,
                                                                       26036, 26099, 38062,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 47522, 0, 3,
                                                                       47018, 37754, 47054,
                                                                       26099, 26162, 38146,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 47630, 0, 3,
                                                                       47054, 37782, 47090,
                                                                       26162, 26225, 38230,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 47738, 0, 3,
                                                                       47090, 37810, 47126,
                                                                       26225, 26288, 38314,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 47846, 0, 3,
                                                                       47126, 37838, 47162,
                                                                       26288, 26351, 38398,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 47954, 0, 3,
                                                                       47198, 37894, 47306,
                                                                       26477, 26603, 38482,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 48170, 0, 3,
                                                                       47306, 37978, 47414,
                                                                       26603, 26729, 38650,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 48386, 0, 3,
                                                                       47414, 38062, 47522,
                                                                       26729, 26855, 38818,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 48602, 0, 3,
                                                                       47522, 38146, 47630,
                                                                       26855, 26981, 38986,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 48818, 0, 3,
                                                                       47630, 38230, 47738,
                                                                       26981, 27107, 39154,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 49034, 0, 3,
                                                                       47738, 38314, 47846,
                                                                       27107, 27233, 39322,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 49250, 0, 3,
                                                                       47954, 38482, 48170,
                                                                       27485, 27695, 39490,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 49610, 0, 3,
                                                                       48170, 38650, 48386,
                                                                       27695, 27905, 39770,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 49970, 0, 3,
                                                                       48386, 38818, 48602,
                                                                       27905, 28115, 40050,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 50330, 0, 3,
                                                                       48602, 38986, 48818,
                                                                       28115, 28325, 40330,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 50690, 0, 3,
                                                                       48818, 39154, 49034,
                                                                       28325, 28535, 40610,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 51050, 0, 3,
                                                                       49250, 39490, 49610,
                                                                       28955, 29270, 40890,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 51590, 0, 3,
                                                                       49610, 39770, 49970,
                                                                       29270, 29585, 41310,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 52130, 0, 3,
                                                                       49970, 40050, 50330,
                                                                       29585, 29900, 41730,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 52670, 0, 3,
                                                                       50330, 40330, 50690,
                                                                       29900, 30215, 42150,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 53210, 0, 3,
                                                                       51050, 40890, 51590,
                                                                       30845, 31286, 42570,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 53966, 0, 3,
                                                                       51590, 41310, 52130,
                                                                       31286, 31727, 43158,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 54722, 0, 3,
                                                                       52130, 41730, 52670,
                                                                       31727, 32168, 43746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 55478, 0, 3,
                                                                       53210, 42570, 53966,
                                                                       33050, 33638, 44334,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 56486, 0, 3,
                                                                       53966, 43158, 54722,
                                                                       33638, 34226, 45118,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 57494, 0, 3,
                                                                       55478, 44334, 56486,
                                                                       35402, 36158, 45902,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 58790, 53210, 756, ncols);

                    simdfunc::contract_primitives(buffer, 59861, 55478, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 61289, 57494, 1296, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 59546, 58790, 21, 1, nmax);

        simdtrf::transform_k_inner(buffer, 60869, 59861, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 62585, 61289, 36, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 63125, 59546, 60869, 15, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 64070, 60869, 62585, 15, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 65330, 63125, 64070, 15, nmax);

        simdtrf::transform_h_inner(buffer, 67220, 65330, 6, 15, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 67220, 165, nmax);
    }

    for (size_t m = 0; m < 825; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
