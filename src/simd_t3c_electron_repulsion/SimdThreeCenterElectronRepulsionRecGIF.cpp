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


#include "SimdThreeCenterElectronRepulsionRecGIF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
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
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_gif_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_gif_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 43015, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 819 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 43015, 22298, 3448, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13}, ncols,
                                                        fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 57, 0, 3, 8, 9,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 63, 0, 3, 9, 10,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 69, 0, 3, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 75, 0, 3, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 81, 0, 3, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 87, 0, 3, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 93, 0, 3, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 99, 0, 3, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 105, 0, 3, 16, 17,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 111, 0, 3, 17, 18,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 117, 0, 3, 18, 19,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 123, 0, 3, 21, 24,
                                                                       57, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 24, 27,
                                                                       63, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 143, 0, 3, 27, 30,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 30, 33,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 33, 36,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 36, 39,
                                                                       87, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 39, 42,
                                                                       93, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 42, 45,
                                                                       99, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 45, 48,
                                                                       105, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 48, 51,
                                                                       111, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 57, 63,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 63, 69,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 69, 75,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 268, 0, 3, 75, 81,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 283, 0, 3, 81, 87,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 87, 93,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 313, 0, 3, 93, 99,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 99,
                                                                       105, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 343, 0, 3, 105,
                                                                       111, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 123,
                                                                       133, 223, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 379, 0, 3, 133,
                                                                       143, 238, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 400, 0, 3, 143,
                                                                       153, 253, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 421, 0, 3, 153,
                                                                       163, 268, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 442, 0, 3, 163,
                                                                       173, 283, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 463, 0, 3, 173,
                                                                       183, 298, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 484, 0, 3, 183,
                                                                       193, 313, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 505, 0, 3, 193,
                                                                       203, 328, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 526, 0, 3, 223,
                                                                       238, 358, 379, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 554, 0, 3, 238,
                                                                       253, 379, 400, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 582, 0, 3, 253,
                                                                       268, 400, 421, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 610, 0, 3, 268,
                                                                       283, 421, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 638, 0, 3, 283,
                                                                       298, 442, 463, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 666, 0, 3, 298,
                                                                       313, 463, 484, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 694, 0, 3, 313,
                                                                       328, 484, 505, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 722, 0, 3, 358,
                                                                       379, 526, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 758, 0, 3, 379,
                                                                       400, 554, 582, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 794, 0, 3, 400,
                                                                       421, 582, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 830, 0, 3, 421,
                                                                       442, 610, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 866, 0, 3, 442,
                                                                       463, 638, 666, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 902, 0, 3, 463,
                                                                       484, 666, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 938, 0, 3, 526,
                                                                       554, 722, 758, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 983, 0, 3, 554,
                                                                       582, 758, 794, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1028, 0, 3, 582,
                                                                       610, 794, 830, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1073, 0, 3, 610,
                                                                       638, 830, 866, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1118, 0, 3, 638,
                                                                       666, 866, 902, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1163, 0, 3, 722,
                                                                       758, 938, 983, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1218, 0, 3, 758,
                                                                       794, 983, 1028, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1273, 0, 3, 794,
                                                                       830, 1028, 1073, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1328, 0, 3, 830,
                                                                       866, 1073, 1118, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1383, 0, 3, 938,
                                                                       983, 1163, 1218, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1449, 0, 3, 983,
                                                                       1028, 1218, 1273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 1515, 0, 3, 1028,
                                                                       1073, 1273, 1328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1581, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1584, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1587, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1590, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1593, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1596, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1599, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1602, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1605, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1608, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1611, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1614, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1617, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1620, 3, 10, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1629, 3, 11, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1638, 3, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1647, 3, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1656, 3, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1665, 3, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1674, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1683, 3, 17, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1692, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1701, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1710, 3, 21, 57,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1728, 3, 24, 63,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1746, 3, 27, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1764, 3, 30, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1782, 3, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1800, 3, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1818, 3, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1836, 3, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1854, 3, 45, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1872, 3, 48, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1890, 3, 51, 117,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1908, 3, 57, 123,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1938, 3, 63, 133,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1968, 3, 69, 143,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1998, 3, 75, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2028, 3, 81, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2058, 3, 87, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2088, 3, 93, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2118, 3, 99, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2148, 3, 105, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2178, 3, 111, 213,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2208, 3, 123, 223,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2253, 3, 133, 238,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2298, 3, 143, 253,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2343, 3, 153, 268,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2388, 3, 163, 283,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2433, 3, 173, 298,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2478, 3, 183, 313,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2523, 3, 193, 328,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2568, 3, 203, 343,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2613, 3, 223, 358,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2676, 3, 238, 379,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2739, 3, 253, 400,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2802, 3, 268, 421,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2865, 3, 283, 442,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2928, 3, 298, 463,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2991, 3, 313, 484,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3054, 3, 328, 505,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3117, 3, 358, 526,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3201, 3, 379, 554,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3285, 3, 400, 582,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3369, 3, 421, 610,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3453, 3, 442, 638,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3537, 3, 463, 666,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3621, 3, 484, 694,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3705, 3, 526, 722,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3813, 3, 554, 758,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3921, 3, 582, 794,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4029, 3, 610, 830,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4137, 3, 638, 866,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4245, 3, 666, 902,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4353, 3, 722, 938,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4488, 3, 758, 983,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4623, 3, 794,
                                                                       1028, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4758, 3, 830,
                                                                       1073, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 4893, 3, 866,
                                                                       1118, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 5028, 3, 938,
                                                                       1163, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 5193, 3, 983,
                                                                       1218, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 5358, 3, 1028,
                                                                       1273, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 5523, 3, 1073,
                                                                       1328, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 5688, 3, 1163,
                                                                       1383, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 5886, 3, 1218,
                                                                       1449, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 6084, 3, 1273,
                                                                       1515, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6282, 3, 8, 9,
                                                                       1587, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6288, 3, 9, 10,
                                                                       1590, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6294, 3, 10, 11,
                                                                       1593, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6300, 3, 11, 12,
                                                                       1596, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6306, 3, 12, 13,
                                                                       1599, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6312, 3, 13, 14,
                                                                       1602, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6318, 3, 14, 15,
                                                                       1605, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6324, 3, 15, 16,
                                                                       1608, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6330, 3, 16, 17,
                                                                       1611, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6336, 3, 17, 18,
                                                                       1614, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6342, 3, 18, 19,
                                                                       1617, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6348, 0, 3, 6282,
                                                                       1587, 6288, 1620, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6366, 0, 3, 6288,
                                                                       1590, 6294, 1629, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6384, 0, 3, 6294,
                                                                       1593, 6300, 1638, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6402, 0, 3, 6300,
                                                                       1596, 6306, 1647, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6420, 0, 3, 6306,
                                                                       1599, 6312, 1656, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6438, 0, 3, 6312,
                                                                       1602, 6318, 1665, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6456, 0, 3, 6318,
                                                                       1605, 6324, 1674, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6474, 0, 3, 6324,
                                                                       1608, 6330, 1683, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6492, 0, 3, 6330,
                                                                       1611, 6336, 1692, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6510, 0, 3, 6336,
                                                                       1614, 6342, 1701, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6528, 0, 3, 6348,
                                                                       1620, 6366, 57, 63, 1746,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6564, 0, 3, 6366,
                                                                       1629, 6384, 63, 69, 1764,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6600, 0, 3, 6384,
                                                                       1638, 6402, 69, 75, 1782,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6636, 0, 3, 6402,
                                                                       1647, 6420, 75, 81, 1800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6672, 0, 3, 6420,
                                                                       1656, 6438, 81, 87, 1818,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6708, 0, 3, 6438,
                                                                       1665, 6456, 87, 93, 1836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6744, 0, 3, 6456,
                                                                       1674, 6474, 93, 99, 1854,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6780, 0, 3, 6474,
                                                                       1683, 6492, 99, 105, 1872,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6816, 0, 3, 6492,
                                                                       1692, 6510, 105, 111,
                                                                       1890, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6852, 0, 3, 6528,
                                                                       1746, 6564, 123, 133,
                                                                       1968, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6912, 0, 3, 6564,
                                                                       1764, 6600, 133, 143,
                                                                       1998, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6972, 0, 3, 6600,
                                                                       1782, 6636, 143, 153,
                                                                       2028, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7032, 0, 3, 6636,
                                                                       1800, 6672, 153, 163,
                                                                       2058, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7092, 0, 3, 6672,
                                                                       1818, 6708, 163, 173,
                                                                       2088, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7152, 0, 3, 6708,
                                                                       1836, 6744, 173, 183,
                                                                       2118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7212, 0, 3, 6744,
                                                                       1854, 6780, 183, 193,
                                                                       2148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7272, 0, 3, 6780,
                                                                       1872, 6816, 193, 203,
                                                                       2178, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7332, 0, 3, 6852,
                                                                       1968, 6912, 223, 238,
                                                                       2298, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7422, 0, 3, 6912,
                                                                       1998, 6972, 238, 253,
                                                                       2343, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7512, 0, 3, 6972,
                                                                       2028, 7032, 253, 268,
                                                                       2388, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7602, 0, 3, 7032,
                                                                       2058, 7092, 268, 283,
                                                                       2433, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7692, 0, 3, 7092,
                                                                       2088, 7152, 283, 298,
                                                                       2478, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7782, 0, 3, 7152,
                                                                       2118, 7212, 298, 313,
                                                                       2523, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7872, 0, 3, 7212,
                                                                       2148, 7272, 313, 328,
                                                                       2568, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7962, 0, 3, 7332,
                                                                       2298, 7422, 358, 379,
                                                                       2739, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8088, 0, 3, 7422,
                                                                       2343, 7512, 379, 400,
                                                                       2802, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8214, 0, 3, 7512,
                                                                       2388, 7602, 400, 421,
                                                                       2865, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8340, 0, 3, 7602,
                                                                       2433, 7692, 421, 442,
                                                                       2928, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8466, 0, 3, 7692,
                                                                       2478, 7782, 442, 463,
                                                                       2991, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 8592, 0, 3, 7782,
                                                                       2523, 7872, 463, 484,
                                                                       3054, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8718, 0, 3, 7962,
                                                                       2739, 8088, 526, 554,
                                                                       3285, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8886, 0, 3, 8088,
                                                                       2802, 8214, 554, 582,
                                                                       3369, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9054, 0, 3, 8214,
                                                                       2865, 8340, 582, 610,
                                                                       3453, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9222, 0, 3, 8340,
                                                                       2928, 8466, 610, 638,
                                                                       3537, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9390, 0, 3, 8466,
                                                                       2991, 8592, 638, 666,
                                                                       3621, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9558, 0, 3, 8718,
                                                                       3285, 8886, 722, 758,
                                                                       3921, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9774, 0, 3, 8886,
                                                                       3369, 9054, 758, 794,
                                                                       4029, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9990, 0, 3, 9054,
                                                                       3453, 9222, 794, 830,
                                                                       4137, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10206, 0, 3, 9222,
                                                                       3537, 9390, 830, 866,
                                                                       4245, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 10422, 0, 3, 9558,
                                                                       3921, 9774, 938, 983,
                                                                       4623, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 10692, 0, 3, 9774,
                                                                       4029, 9990, 983, 1028,
                                                                       4758, ncols, gamma, p,
                                                                       q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 10962, 0, 3, 9990,
                                                                       4137, 10206, 1028, 1073,
                                                                       4893, ncols, gamma, p,
                                                                       q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 11232, 0, 3,
                                                                       10422, 4623, 10692, 1163,
                                                                       1218, 5358, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 11562, 0, 3,
                                                                       10692, 4758, 10962, 1218,
                                                                       1273, 5523, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 11892, 0, 3,
                                                                       11232, 5358, 11562, 1383,
                                                                       1449, 6084, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12288, 3, 1581,
                                                                       1584, 6282, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12298, 3, 1584,
                                                                       1587, 6288, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12308, 3, 1587,
                                                                       1590, 6294, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12318, 3, 1590,
                                                                       1593, 6300, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12328, 3, 1593,
                                                                       1596, 6306, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12338, 3, 1596,
                                                                       1599, 6312, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12348, 3, 1599,
                                                                       1602, 6318, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12358, 3, 1602,
                                                                       1605, 6324, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12368, 3, 1605,
                                                                       1608, 6330, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12378, 3, 1608,
                                                                       1611, 6336, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12388, 3, 1611,
                                                                       1614, 6342, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12398, 0, 3,
                                                                       12288, 6282, 12298, 6348,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12428, 0, 3,
                                                                       12298, 6288, 12308, 6366,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12458, 0, 3,
                                                                       12308, 6294, 12318, 6384,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12488, 0, 3,
                                                                       12318, 6300, 12328, 6402,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12518, 0, 3,
                                                                       12328, 6306, 12338, 6420,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12548, 0, 3,
                                                                       12338, 6312, 12348, 6438,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12578, 0, 3,
                                                                       12348, 6318, 12358, 6456,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12608, 0, 3,
                                                                       12358, 6324, 12368, 6474,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12638, 0, 3,
                                                                       12368, 6330, 12378, 6492,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12668, 0, 3,
                                                                       12378, 6336, 12388, 6510,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12698, 0, 3,
                                                                       12398, 6348, 12428, 1710,
                                                                       1728, 6528, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12758, 0, 3,
                                                                       12428, 6366, 12458, 1728,
                                                                       1746, 6564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12818, 0, 3,
                                                                       12458, 6384, 12488, 1746,
                                                                       1764, 6600, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12878, 0, 3,
                                                                       12488, 6402, 12518, 1764,
                                                                       1782, 6636, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12938, 0, 3,
                                                                       12518, 6420, 12548, 1782,
                                                                       1800, 6672, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12998, 0, 3,
                                                                       12548, 6438, 12578, 1800,
                                                                       1818, 6708, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13058, 0, 3,
                                                                       12578, 6456, 12608, 1818,
                                                                       1836, 6744, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13118, 0, 3,
                                                                       12608, 6474, 12638, 1836,
                                                                       1854, 6780, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13178, 0, 3,
                                                                       12638, 6492, 12668, 1854,
                                                                       1872, 6816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13238, 0, 3,
                                                                       12698, 6528, 12758, 1908,
                                                                       1938, 6852, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13338, 0, 3,
                                                                       12758, 6564, 12818, 1938,
                                                                       1968, 6912, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13438, 0, 3,
                                                                       12818, 6600, 12878, 1968,
                                                                       1998, 6972, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13538, 0, 3,
                                                                       12878, 6636, 12938, 1998,
                                                                       2028, 7032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13638, 0, 3,
                                                                       12938, 6672, 12998, 2028,
                                                                       2058, 7092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13738, 0, 3,
                                                                       12998, 6708, 13058, 2058,
                                                                       2088, 7152, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13838, 0, 3,
                                                                       13058, 6744, 13118, 2088,
                                                                       2118, 7212, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13938, 0, 3,
                                                                       13118, 6780, 13178, 2118,
                                                                       2148, 7272, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14038, 0, 3,
                                                                       13238, 6852, 13338, 2208,
                                                                       2253, 7332, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14188, 0, 3,
                                                                       13338, 6912, 13438, 2253,
                                                                       2298, 7422, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14338, 0, 3,
                                                                       13438, 6972, 13538, 2298,
                                                                       2343, 7512, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14488, 0, 3,
                                                                       13538, 7032, 13638, 2343,
                                                                       2388, 7602, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14638, 0, 3,
                                                                       13638, 7092, 13738, 2388,
                                                                       2433, 7692, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14788, 0, 3,
                                                                       13738, 7152, 13838, 2433,
                                                                       2478, 7782, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14938, 0, 3,
                                                                       13838, 7212, 13938, 2478,
                                                                       2523, 7872, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15088, 0, 3,
                                                                       14038, 7332, 14188, 2613,
                                                                       2676, 7962, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15298, 0, 3,
                                                                       14188, 7422, 14338, 2676,
                                                                       2739, 8088, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15508, 0, 3,
                                                                       14338, 7512, 14488, 2739,
                                                                       2802, 8214, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15718, 0, 3,
                                                                       14488, 7602, 14638, 2802,
                                                                       2865, 8340, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15928, 0, 3,
                                                                       14638, 7692, 14788, 2865,
                                                                       2928, 8466, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 16138, 0, 3,
                                                                       14788, 7782, 14938, 2928,
                                                                       2991, 8592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16348, 0, 3,
                                                                       15088, 7962, 15298, 3117,
                                                                       3201, 8718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16628, 0, 3,
                                                                       15298, 8088, 15508, 3201,
                                                                       3285, 8886, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16908, 0, 3,
                                                                       15508, 8214, 15718, 3285,
                                                                       3369, 9054, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 17188, 0, 3,
                                                                       15718, 8340, 15928, 3369,
                                                                       3453, 9222, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 17468, 0, 3,
                                                                       15928, 8466, 16138, 3453,
                                                                       3537, 9390, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 17748, 0, 3,
                                                                       16348, 8718, 16628, 3705,
                                                                       3813, 9558, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 18108, 0, 3,
                                                                       16628, 8886, 16908, 3813,
                                                                       3921, 9774, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 18468, 0, 3,
                                                                       16908, 9054, 17188, 3921,
                                                                       4029, 9990, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 18828, 0, 3,
                                                                       17188, 9222, 17468, 4029,
                                                                       4137, 10206, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 19188, 0, 3,
                                                                       17748, 9558, 18108, 4353,
                                                                       4488, 10422, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 19638, 0, 3,
                                                                       18108, 9774, 18468, 4488,
                                                                       4623, 10692, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 20088, 0, 3,
                                                                       18468, 9990, 18828, 4623,
                                                                       4758, 10962, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 20538, 0, 3,
                                                                       19188, 10422, 19638, 5028,
                                                                       5193, 11232, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 21088, 0, 3,
                                                                       19638, 10692, 20088, 5193,
                                                                       5358, 11562, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 21638, 0, 3,
                                                                       20538, 11232, 21088, 5688,
                                                                       5886, 11892, ncols, gamma,
                                                                       p, q);

                    simdfunc::contract_primitives(buffer, 22298, 16348, 280, ncols);

                    simdfunc::contract_primitives(buffer, 22774, 17748, 360, ncols);

                    simdfunc::contract_primitives(buffer, 23386, 19188, 450, ncols);

                    simdfunc::contract_primitives(buffer, 24151, 20538, 550, ncols);

                    simdfunc::contract_primitives(buffer, 25086, 21638, 660, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 22578, 22298, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 23134, 22774, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 23836, 23386, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 24701, 24151, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 25746, 25086, 66, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 26208, 22578, 23134, 7, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 26796, 23134, 23836, 7, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 27552, 23836, 24701, 7, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 28497, 24701, 25746, 7, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 29652, 26208, 26796, 7, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 30828, 26796, 27552, 7, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 32340, 27552, 28497, 7, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 34230, 29652, 30828, 7, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 36190, 30828, 32340, 7, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 38710, 34230, 36190, 7, nmax);

        simdtrf::transform_i_inner(buffer, 41650, 38710, 15, 7, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 41650, 91, nmax);
    }

    for (size_t m = 0; m < 819; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
