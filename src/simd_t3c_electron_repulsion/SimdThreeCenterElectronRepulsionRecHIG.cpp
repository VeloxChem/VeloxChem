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


#include "SimdThreeCenterElectronRepulsionRecHIG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hig_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hig_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 104448, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1287 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 104448, 55791, 6690, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 15,
                                                             ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 69, 0, 3, 8, 9,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 75, 0, 3, 9, 10,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 81, 0, 3, 10, 11,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 87, 0, 3, 11, 12,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 93, 0, 3, 12, 13,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 99, 0, 3, 13, 14,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 105, 0, 3, 14, 15,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 111, 0, 3, 15, 16,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 117, 0, 3, 16, 17,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 123, 0, 3, 17, 18,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 129, 0, 3, 18, 19,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 135, 0, 3, 19, 20,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 141, 0, 3, 20, 21,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 147, 0, 3, 21, 22,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 153, 0, 3, 24, 27,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 163, 0, 3, 27, 30,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 30, 33,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 183, 0, 3, 33, 36,
                                                                       87, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 193, 0, 3, 36, 39,
                                                                       93, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 39, 42,
                                                                       99, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 42, 45,
                                                                       105, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 223, 0, 3, 45, 48,
                                                                       111, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 233, 0, 3, 48, 51,
                                                                       117, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 51, 54,
                                                                       123, 129, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 253, 0, 3, 54, 57,
                                                                       129, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 263, 0, 3, 57, 60,
                                                                       135, 141, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 60, 63,
                                                                       141, 147, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 283, 0, 3, 69, 75,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 298, 0, 3, 75, 81,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 313, 0, 3, 81, 87,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 328, 0, 3, 87, 93,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 343, 0, 3, 93, 99,
                                                                       193, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 358, 0, 3, 99,
                                                                       105, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 373, 0, 3, 105,
                                                                       111, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 388, 0, 3, 111,
                                                                       117, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 403, 0, 3, 117,
                                                                       123, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 418, 0, 3, 123,
                                                                       129, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 433, 0, 3, 129,
                                                                       135, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 448, 0, 3, 135,
                                                                       141, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 463, 0, 3, 153,
                                                                       163, 283, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 484, 0, 3, 163,
                                                                       173, 298, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 505, 0, 3, 173,
                                                                       183, 313, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 526, 0, 3, 183,
                                                                       193, 328, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 547, 0, 3, 193,
                                                                       203, 343, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 568, 0, 3, 203,
                                                                       213, 358, 373, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 589, 0, 3, 213,
                                                                       223, 373, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 610, 0, 3, 223,
                                                                       233, 388, 403, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 631, 0, 3, 233,
                                                                       243, 403, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 652, 0, 3, 243,
                                                                       253, 418, 433, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 673, 0, 3, 253,
                                                                       263, 433, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 694, 0, 3, 283,
                                                                       298, 463, 484, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 722, 0, 3, 298,
                                                                       313, 484, 505, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 750, 0, 3, 313,
                                                                       328, 505, 526, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 778, 0, 3, 328,
                                                                       343, 526, 547, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 806, 0, 3, 343,
                                                                       358, 547, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 834, 0, 3, 358,
                                                                       373, 568, 589, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 862, 0, 3, 373,
                                                                       388, 589, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 890, 0, 3, 388,
                                                                       403, 610, 631, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 918, 0, 3, 403,
                                                                       418, 631, 652, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 946, 0, 3, 418,
                                                                       433, 652, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 974, 0, 3, 463,
                                                                       484, 694, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1010, 0, 3, 484,
                                                                       505, 722, 750, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1046, 0, 3, 505,
                                                                       526, 750, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1082, 0, 3, 526,
                                                                       547, 778, 806, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1118, 0, 3, 547,
                                                                       568, 806, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1154, 0, 3, 568,
                                                                       589, 834, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1190, 0, 3, 589,
                                                                       610, 862, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1226, 0, 3, 610,
                                                                       631, 890, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 1262, 0, 3, 631,
                                                                       652, 918, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1298, 0, 3, 694,
                                                                       722, 974, 1010, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1343, 0, 3, 722,
                                                                       750, 1010, 1046, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1388, 0, 3, 750,
                                                                       778, 1046, 1082, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1433, 0, 3, 778,
                                                                       806, 1082, 1118, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 806,
                                                                       834, 1118, 1154, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1523, 0, 3, 834,
                                                                       862, 1154, 1190, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1568, 0, 3, 862,
                                                                       890, 1190, 1226, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 1613, 0, 3, 890,
                                                                       918, 1226, 1262, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1658, 0, 3, 974,
                                                                       1010, 1298, 1343, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1713, 0, 3, 1010,
                                                                       1046, 1343, 1388, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1768, 0, 3, 1046,
                                                                       1082, 1388, 1433, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1823, 0, 3, 1082,
                                                                       1118, 1433, 1478, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1878, 0, 3, 1118,
                                                                       1154, 1478, 1523, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1933, 0, 3, 1154,
                                                                       1190, 1523, 1568, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1190,
                                                                       1226, 1568, 1613, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2043, 0, 3, 1298,
                                                                       1343, 1658, 1713, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2109, 0, 3, 1343,
                                                                       1388, 1713, 1768, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2175, 0, 3, 1388,
                                                                       1433, 1768, 1823, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2241, 0, 3, 1433,
                                                                       1478, 1823, 1878, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2307, 0, 3, 1478,
                                                                       1523, 1878, 1933, ncols,
                                                                       gamma, p, q);

                    compute_prim_sns_three_center_electron_repulsion_0(buffer, 2373, 0, 3, 1523,
                                                                       1568, 1933, 1988, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2439, 0, 3, 1658,
                                                                       1713, 2043, 2109, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2517, 0, 3, 1713,
                                                                       1768, 2109, 2175, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2595, 0, 3, 1768,
                                                                       1823, 2175, 2241, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2673, 0, 3, 1823,
                                                                       1878, 2241, 2307, ncols,
                                                                       gamma, p, q);

                    compute_prim_sos_three_center_electron_repulsion_0(buffer, 2751, 0, 3, 1878,
                                                                       1933, 2307, 2373, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2829, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2832, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2835, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2838, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2841, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2844, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2847, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2850, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2853, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2856, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2859, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2862, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2865, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2868, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2871, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2880, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2889, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2898, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2907, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2916, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2925, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2934, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2943, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2952, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2961, 3, 20, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2970, 3, 21, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2979, 3, 22, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2988, 3, 30, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3006, 3, 33, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3024, 3, 36, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3042, 3, 39, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3060, 3, 42, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3078, 3, 45, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3096, 3, 48, 117,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3114, 3, 51, 123,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3132, 3, 54, 129,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3150, 3, 57, 135,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3168, 3, 60, 141,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3186, 3, 63, 147,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3204, 3, 81, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3234, 3, 87, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3264, 3, 93, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3294, 3, 99, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3324, 3, 105, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3354, 3, 111, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3384, 3, 117, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3414, 3, 123, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3444, 3, 129, 253,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3474, 3, 135, 263,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3504, 3, 141, 273,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3534, 3, 173, 313,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3579, 3, 183, 328,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3624, 3, 193, 343,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3669, 3, 203, 358,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3714, 3, 213, 373,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3759, 3, 223, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3804, 3, 233, 403,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3849, 3, 243, 418,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3894, 3, 253, 433,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3939, 3, 263, 448,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3984, 3, 313, 505,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4047, 3, 328, 526,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4110, 3, 343, 547,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4173, 3, 358, 568,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4236, 3, 373, 589,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4299, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4362, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4425, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4488, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4551, 3, 505, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4635, 3, 526, 778,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4719, 3, 547, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4803, 3, 568, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4887, 3, 589, 862,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4971, 3, 610, 890,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5055, 3, 631, 918,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5139, 3, 652, 946,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5223, 3, 750,
                                                                       1046, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5331, 3, 778,
                                                                       1082, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5439, 3, 806,
                                                                       1118, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5547, 3, 834,
                                                                       1154, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5655, 3, 862,
                                                                       1190, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5763, 3, 890,
                                                                       1226, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5871, 3, 918,
                                                                       1262, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5979, 3, 1046,
                                                                       1388, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6114, 3, 1082,
                                                                       1433, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6249, 3, 1118,
                                                                       1478, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6384, 3, 1154,
                                                                       1523, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6519, 3, 1190,
                                                                       1568, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6654, 3, 1226,
                                                                       1613, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6789, 3, 1388,
                                                                       1768, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6954, 3, 1433,
                                                                       1823, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7119, 3, 1478,
                                                                       1878, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7284, 3, 1523,
                                                                       1933, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7449, 3, 1568,
                                                                       1988, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7614, 3, 1768,
                                                                       2175, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 7812, 3, 1823,
                                                                       2241, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8010, 3, 1878,
                                                                       2307, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8208, 3, 1933,
                                                                       2373, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 8406, 3, 2175,
                                                                       2595, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 8640, 3, 2241,
                                                                       2673, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 8874, 3, 2307,
                                                                       2751, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9108, 3, 8, 9,
                                                                       2829, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9114, 3, 9, 10,
                                                                       2832, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9120, 3, 10, 11,
                                                                       2835, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9126, 3, 11, 12,
                                                                       2838, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9132, 3, 12, 13,
                                                                       2841, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9138, 3, 13, 14,
                                                                       2844, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9144, 3, 14, 15,
                                                                       2847, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9150, 3, 15, 16,
                                                                       2850, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9156, 3, 16, 17,
                                                                       2853, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9162, 3, 17, 18,
                                                                       2856, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9168, 3, 18, 19,
                                                                       2859, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9174, 3, 19, 20,
                                                                       2862, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9180, 3, 20, 21,
                                                                       2865, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 9186, 3, 21, 22,
                                                                       2868, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9192, 0, 3, 9108,
                                                                       2829, 9114, 2871, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9210, 0, 3, 9114,
                                                                       2832, 9120, 2880, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9228, 0, 3, 9120,
                                                                       2835, 9126, 2889, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9246, 0, 3, 9126,
                                                                       2838, 9132, 2898, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9264, 0, 3, 9132,
                                                                       2841, 9138, 2907, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9282, 0, 3, 9138,
                                                                       2844, 9144, 2916, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9300, 0, 3, 9144,
                                                                       2847, 9150, 2925, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9318, 0, 3, 9150,
                                                                       2850, 9156, 2934, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9336, 0, 3, 9156,
                                                                       2853, 9162, 2943, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9354, 0, 3, 9162,
                                                                       2856, 9168, 2952, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9372, 0, 3, 9168,
                                                                       2859, 9174, 2961, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9390, 0, 3, 9174,
                                                                       2862, 9180, 2970, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 9408, 0, 3, 9180,
                                                                       2865, 9186, 2979, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9426, 0, 3, 9192,
                                                                       2871, 9210, 69, 75, 2988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9462, 0, 3, 9210,
                                                                       2880, 9228, 75, 81, 3006,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9498, 0, 3, 9228,
                                                                       2889, 9246, 81, 87, 3024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9534, 0, 3, 9246,
                                                                       2898, 9264, 87, 93, 3042,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9570, 0, 3, 9264,
                                                                       2907, 9282, 93, 99, 3060,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9606, 0, 3, 9282,
                                                                       2916, 9300, 99, 105, 3078,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9642, 0, 3, 9300,
                                                                       2925, 9318, 105, 111,
                                                                       3096, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9678, 0, 3, 9318,
                                                                       2934, 9336, 111, 117,
                                                                       3114, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9714, 0, 3, 9336,
                                                                       2943, 9354, 117, 123,
                                                                       3132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9750, 0, 3, 9354,
                                                                       2952, 9372, 123, 129,
                                                                       3150, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9786, 0, 3, 9372,
                                                                       2961, 9390, 129, 135,
                                                                       3168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 9822, 0, 3, 9390,
                                                                       2970, 9408, 135, 141,
                                                                       3186, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9858, 0, 3, 9426,
                                                                       2988, 9462, 153, 163,
                                                                       3204, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9918, 0, 3, 9462,
                                                                       3006, 9498, 163, 173,
                                                                       3234, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9978, 0, 3, 9498,
                                                                       3024, 9534, 173, 183,
                                                                       3264, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10038, 0, 3, 9534,
                                                                       3042, 9570, 183, 193,
                                                                       3294, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10098, 0, 3, 9570,
                                                                       3060, 9606, 193, 203,
                                                                       3324, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10158, 0, 3, 9606,
                                                                       3078, 9642, 203, 213,
                                                                       3354, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10218, 0, 3, 9642,
                                                                       3096, 9678, 213, 223,
                                                                       3384, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10278, 0, 3, 9678,
                                                                       3114, 9714, 223, 233,
                                                                       3414, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10338, 0, 3, 9714,
                                                                       3132, 9750, 233, 243,
                                                                       3444, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10398, 0, 3, 9750,
                                                                       3150, 9786, 243, 253,
                                                                       3474, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 10458, 0, 3, 9786,
                                                                       3168, 9822, 253, 263,
                                                                       3504, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10518, 0, 3, 9858,
                                                                       3204, 9918, 283, 298,
                                                                       3534, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10608, 0, 3, 9918,
                                                                       3234, 9978, 298, 313,
                                                                       3579, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10698, 0, 3, 9978,
                                                                       3264, 10038, 313, 328,
                                                                       3624, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10788, 0, 3,
                                                                       10038, 3294, 10098, 328,
                                                                       343, 3669, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10878, 0, 3,
                                                                       10098, 3324, 10158, 343,
                                                                       358, 3714, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10968, 0, 3,
                                                                       10158, 3354, 10218, 358,
                                                                       373, 3759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11058, 0, 3,
                                                                       10218, 3384, 10278, 373,
                                                                       388, 3804, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11148, 0, 3,
                                                                       10278, 3414, 10338, 388,
                                                                       403, 3849, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11238, 0, 3,
                                                                       10338, 3444, 10398, 403,
                                                                       418, 3894, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 11328, 0, 3,
                                                                       10398, 3474, 10458, 418,
                                                                       433, 3939, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11418, 0, 3,
                                                                       10518, 3534, 10608, 463,
                                                                       484, 3984, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11544, 0, 3,
                                                                       10608, 3579, 10698, 484,
                                                                       505, 4047, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11670, 0, 3,
                                                                       10698, 3624, 10788, 505,
                                                                       526, 4110, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11796, 0, 3,
                                                                       10788, 3669, 10878, 526,
                                                                       547, 4173, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11922, 0, 3,
                                                                       10878, 3714, 10968, 547,
                                                                       568, 4236, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12048, 0, 3,
                                                                       10968, 3759, 11058, 568,
                                                                       589, 4299, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12174, 0, 3,
                                                                       11058, 3804, 11148, 589,
                                                                       610, 4362, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12300, 0, 3,
                                                                       11148, 3849, 11238, 610,
                                                                       631, 4425, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 12426, 0, 3,
                                                                       11238, 3894, 11328, 631,
                                                                       652, 4488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12552, 0, 3,
                                                                       11418, 3984, 11544, 694,
                                                                       722, 4551, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12720, 0, 3,
                                                                       11544, 4047, 11670, 722,
                                                                       750, 4635, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12888, 0, 3,
                                                                       11670, 4110, 11796, 750,
                                                                       778, 4719, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13056, 0, 3,
                                                                       11796, 4173, 11922, 778,
                                                                       806, 4803, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13224, 0, 3,
                                                                       11922, 4236, 12048, 806,
                                                                       834, 4887, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13392, 0, 3,
                                                                       12048, 4299, 12174, 834,
                                                                       862, 4971, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13560, 0, 3,
                                                                       12174, 4362, 12300, 862,
                                                                       890, 5055, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 13728, 0, 3,
                                                                       12300, 4425, 12426, 890,
                                                                       918, 5139, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13896, 0, 3,
                                                                       12552, 4551, 12720, 974,
                                                                       1010, 5223, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14112, 0, 3,
                                                                       12720, 4635, 12888, 1010,
                                                                       1046, 5331, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14328, 0, 3,
                                                                       12888, 4719, 13056, 1046,
                                                                       1082, 5439, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14544, 0, 3,
                                                                       13056, 4803, 13224, 1082,
                                                                       1118, 5547, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14760, 0, 3,
                                                                       13224, 4887, 13392, 1118,
                                                                       1154, 5655, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14976, 0, 3,
                                                                       13392, 4971, 13560, 1154,
                                                                       1190, 5763, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 15192, 0, 3,
                                                                       13560, 5055, 13728, 1190,
                                                                       1226, 5871, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15408, 0, 3,
                                                                       13896, 5223, 14112, 1298,
                                                                       1343, 5979, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15678, 0, 3,
                                                                       14112, 5331, 14328, 1343,
                                                                       1388, 6114, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15948, 0, 3,
                                                                       14328, 5439, 14544, 1388,
                                                                       1433, 6249, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16218, 0, 3,
                                                                       14544, 5547, 14760, 1433,
                                                                       1478, 6384, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16488, 0, 3,
                                                                       14760, 5655, 14976, 1478,
                                                                       1523, 6519, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 16758, 0, 3,
                                                                       14976, 5763, 15192, 1523,
                                                                       1568, 6654, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 17028, 0, 3,
                                                                       15408, 5979, 15678, 1658,
                                                                       1713, 6789, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 17358, 0, 3,
                                                                       15678, 6114, 15948, 1713,
                                                                       1768, 6954, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 17688, 0, 3,
                                                                       15948, 6249, 16218, 1768,
                                                                       1823, 7119, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 18018, 0, 3,
                                                                       16218, 6384, 16488, 1823,
                                                                       1878, 7284, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 18348, 0, 3,
                                                                       16488, 6519, 16758, 1878,
                                                                       1933, 7449, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 18678, 0, 3,
                                                                       17028, 6789, 17358, 2043,
                                                                       2109, 7614, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 19074, 0, 3,
                                                                       17358, 6954, 17688, 2109,
                                                                       2175, 7812, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 19470, 0, 3,
                                                                       17688, 7119, 18018, 2175,
                                                                       2241, 8010, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 19866, 0, 3,
                                                                       18018, 7284, 18348, 2241,
                                                                       2307, 8208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 20262, 0, 3,
                                                                       18678, 7614, 19074, 2439,
                                                                       2517, 8406, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 20730, 0, 3,
                                                                       19074, 7812, 19470, 2517,
                                                                       2595, 8640, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 21198, 0, 3,
                                                                       19470, 8010, 19866, 2595,
                                                                       2673, 8874, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21666, 3, 2829,
                                                                       2832, 9120, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21676, 3, 2832,
                                                                       2835, 9126, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21686, 3, 2835,
                                                                       2838, 9132, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21696, 3, 2838,
                                                                       2841, 9138, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21706, 3, 2841,
                                                                       2844, 9144, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21716, 3, 2844,
                                                                       2847, 9150, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21726, 3, 2847,
                                                                       2850, 9156, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21736, 3, 2850,
                                                                       2853, 9162, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21746, 3, 2853,
                                                                       2856, 9168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21756, 3, 2856,
                                                                       2859, 9174, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21766, 3, 2859,
                                                                       2862, 9180, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21776, 3, 2862,
                                                                       2865, 9186, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21786, 0, 3,
                                                                       21666, 9120, 21676, 9228,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21816, 0, 3,
                                                                       21676, 9126, 21686, 9246,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21846, 0, 3,
                                                                       21686, 9132, 21696, 9264,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21876, 0, 3,
                                                                       21696, 9138, 21706, 9282,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21906, 0, 3,
                                                                       21706, 9144, 21716, 9300,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21936, 0, 3,
                                                                       21716, 9150, 21726, 9318,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21966, 0, 3,
                                                                       21726, 9156, 21736, 9336,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 21996, 0, 3,
                                                                       21736, 9162, 21746, 9354,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22026, 0, 3,
                                                                       21746, 9168, 21756, 9372,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22056, 0, 3,
                                                                       21756, 9174, 21766, 9390,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 22086, 0, 3,
                                                                       21766, 9180, 21776, 9408,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22116, 0, 3,
                                                                       21786, 9228, 21816, 2988,
                                                                       3006, 9498, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22176, 0, 3,
                                                                       21816, 9246, 21846, 3006,
                                                                       3024, 9534, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22236, 0, 3,
                                                                       21846, 9264, 21876, 3024,
                                                                       3042, 9570, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22296, 0, 3,
                                                                       21876, 9282, 21906, 3042,
                                                                       3060, 9606, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22356, 0, 3,
                                                                       21906, 9300, 21936, 3060,
                                                                       3078, 9642, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22416, 0, 3,
                                                                       21936, 9318, 21966, 3078,
                                                                       3096, 9678, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22476, 0, 3,
                                                                       21966, 9336, 21996, 3096,
                                                                       3114, 9714, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22536, 0, 3,
                                                                       21996, 9354, 22026, 3114,
                                                                       3132, 9750, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22596, 0, 3,
                                                                       22026, 9372, 22056, 3132,
                                                                       3150, 9786, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 22656, 0, 3,
                                                                       22056, 9390, 22086, 3150,
                                                                       3168, 9822, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22716, 0, 3,
                                                                       22116, 9498, 22176, 3204,
                                                                       3234, 9978, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22816, 0, 3,
                                                                       22176, 9534, 22236, 3234,
                                                                       3264, 10038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 22916, 0, 3,
                                                                       22236, 9570, 22296, 3264,
                                                                       3294, 10098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23016, 0, 3,
                                                                       22296, 9606, 22356, 3294,
                                                                       3324, 10158, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23116, 0, 3,
                                                                       22356, 9642, 22416, 3324,
                                                                       3354, 10218, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23216, 0, 3,
                                                                       22416, 9678, 22476, 3354,
                                                                       3384, 10278, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23316, 0, 3,
                                                                       22476, 9714, 22536, 3384,
                                                                       3414, 10338, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23416, 0, 3,
                                                                       22536, 9750, 22596, 3414,
                                                                       3444, 10398, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 23516, 0, 3,
                                                                       22596, 9786, 22656, 3444,
                                                                       3474, 10458, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23616, 0, 3,
                                                                       22716, 9978, 22816, 3534,
                                                                       3579, 10698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23766, 0, 3,
                                                                       22816, 10038, 22916, 3579,
                                                                       3624, 10788, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 23916, 0, 3,
                                                                       22916, 10098, 23016, 3624,
                                                                       3669, 10878, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24066, 0, 3,
                                                                       23016, 10158, 23116, 3669,
                                                                       3714, 10968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24216, 0, 3,
                                                                       23116, 10218, 23216, 3714,
                                                                       3759, 11058, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24366, 0, 3,
                                                                       23216, 10278, 23316, 3759,
                                                                       3804, 11148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24516, 0, 3,
                                                                       23316, 10338, 23416, 3804,
                                                                       3849, 11238, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 24666, 0, 3,
                                                                       23416, 10398, 23516, 3849,
                                                                       3894, 11328, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 24816, 0, 3,
                                                                       23616, 10698, 23766, 3984,
                                                                       4047, 11670, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25026, 0, 3,
                                                                       23766, 10788, 23916, 4047,
                                                                       4110, 11796, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25236, 0, 3,
                                                                       23916, 10878, 24066, 4110,
                                                                       4173, 11922, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25446, 0, 3,
                                                                       24066, 10968, 24216, 4173,
                                                                       4236, 12048, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25656, 0, 3,
                                                                       24216, 11058, 24366, 4236,
                                                                       4299, 12174, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 25866, 0, 3,
                                                                       24366, 11148, 24516, 4299,
                                                                       4362, 12300, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 26076, 0, 3,
                                                                       24516, 11238, 24666, 4362,
                                                                       4425, 12426, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 26286, 0, 3,
                                                                       24816, 11670, 25026, 4551,
                                                                       4635, 12888, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 26566, 0, 3,
                                                                       25026, 11796, 25236, 4635,
                                                                       4719, 13056, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 26846, 0, 3,
                                                                       25236, 11922, 25446, 4719,
                                                                       4803, 13224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27126, 0, 3,
                                                                       25446, 12048, 25656, 4803,
                                                                       4887, 13392, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27406, 0, 3,
                                                                       25656, 12174, 25866, 4887,
                                                                       4971, 13560, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 27686, 0, 3,
                                                                       25866, 12300, 26076, 4971,
                                                                       5055, 13728, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 27966, 0, 3,
                                                                       26286, 12888, 26566, 5223,
                                                                       5331, 14328, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 28326, 0, 3,
                                                                       26566, 13056, 26846, 5331,
                                                                       5439, 14544, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 28686, 0, 3,
                                                                       26846, 13224, 27126, 5439,
                                                                       5547, 14760, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29046, 0, 3,
                                                                       27126, 13392, 27406, 5547,
                                                                       5655, 14976, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 29406, 0, 3,
                                                                       27406, 13560, 27686, 5655,
                                                                       5763, 15192, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 29766, 0, 3,
                                                                       27966, 14328, 28326, 5979,
                                                                       6114, 15948, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 30216, 0, 3,
                                                                       28326, 14544, 28686, 6114,
                                                                       6249, 16218, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 30666, 0, 3,
                                                                       28686, 14760, 29046, 6249,
                                                                       6384, 16488, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 31116, 0, 3,
                                                                       29046, 14976, 29406, 6384,
                                                                       6519, 16758, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 31566, 0, 3,
                                                                       29766, 15948, 30216, 6789,
                                                                       6954, 17688, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 32116, 0, 3,
                                                                       30216, 16218, 30666, 6954,
                                                                       7119, 18018, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 32666, 0, 3,
                                                                       30666, 16488, 31116, 7119,
                                                                       7284, 18348, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 33216, 0, 3,
                                                                       31566, 17688, 32116, 7614,
                                                                       7812, 19470, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 33876, 0, 3,
                                                                       32116, 18018, 32666, 7812,
                                                                       8010, 19866, ncols, gamma,
                                                                       p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 34536, 0, 3,
                                                                       33216, 19470, 33876, 8406,
                                                                       8640, 21198, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35316, 3, 9108,
                                                                       9114, 21666, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35331, 3, 9114,
                                                                       9120, 21676, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35346, 3, 9120,
                                                                       9126, 21686, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35361, 3, 9126,
                                                                       9132, 21696, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35376, 3, 9132,
                                                                       9138, 21706, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35391, 3, 9138,
                                                                       9144, 21716, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35406, 3, 9144,
                                                                       9150, 21726, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35421, 3, 9150,
                                                                       9156, 21736, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35436, 3, 9156,
                                                                       9162, 21746, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35451, 3, 9162,
                                                                       9168, 21756, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35466, 3, 9168,
                                                                       9174, 21766, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 35481, 3, 9174,
                                                                       9180, 21776, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35496, 0, 3,
                                                                       35316, 21666, 35331, 9192,
                                                                       9210, 21786, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35541, 0, 3,
                                                                       35331, 21676, 35346, 9210,
                                                                       9228, 21816, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35586, 0, 3,
                                                                       35346, 21686, 35361, 9228,
                                                                       9246, 21846, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35631, 0, 3,
                                                                       35361, 21696, 35376, 9246,
                                                                       9264, 21876, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35676, 0, 3,
                                                                       35376, 21706, 35391, 9264,
                                                                       9282, 21906, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35721, 0, 3,
                                                                       35391, 21716, 35406, 9282,
                                                                       9300, 21936, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35766, 0, 3,
                                                                       35406, 21726, 35421, 9300,
                                                                       9318, 21966, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35811, 0, 3,
                                                                       35421, 21736, 35436, 9318,
                                                                       9336, 21996, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35856, 0, 3,
                                                                       35436, 21746, 35451, 9336,
                                                                       9354, 22026, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35901, 0, 3,
                                                                       35451, 21756, 35466, 9354,
                                                                       9372, 22056, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 35946, 0, 3,
                                                                       35466, 21766, 35481, 9372,
                                                                       9390, 22086, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35991, 0, 3,
                                                                       35496, 21786, 35541, 9426,
                                                                       9462, 22116, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36081, 0, 3,
                                                                       35541, 21816, 35586, 9462,
                                                                       9498, 22176, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36171, 0, 3,
                                                                       35586, 21846, 35631, 9498,
                                                                       9534, 22236, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36261, 0, 3,
                                                                       35631, 21876, 35676, 9534,
                                                                       9570, 22296, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36351, 0, 3,
                                                                       35676, 21906, 35721, 9570,
                                                                       9606, 22356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36441, 0, 3,
                                                                       35721, 21936, 35766, 9606,
                                                                       9642, 22416, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36531, 0, 3,
                                                                       35766, 21966, 35811, 9642,
                                                                       9678, 22476, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36621, 0, 3,
                                                                       35811, 21996, 35856, 9678,
                                                                       9714, 22536, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36711, 0, 3,
                                                                       35856, 22026, 35901, 9714,
                                                                       9750, 22596, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 36801, 0, 3,
                                                                       35901, 22056, 35946, 9750,
                                                                       9786, 22656, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36891, 0, 3,
                                                                       35991, 22116, 36081, 9858,
                                                                       9918, 22716, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37041, 0, 3,
                                                                       36081, 22176, 36171, 9918,
                                                                       9978, 22816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37191, 0, 3,
                                                                       36171, 22236, 36261, 9978,
                                                                       10038, 22916, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37341, 0, 3,
                                                                       36261, 22296, 36351,
                                                                       10038, 10098, 23016,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37491, 0, 3,
                                                                       36351, 22356, 36441,
                                                                       10098, 10158, 23116,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37641, 0, 3,
                                                                       36441, 22416, 36531,
                                                                       10158, 10218, 23216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37791, 0, 3,
                                                                       36531, 22476, 36621,
                                                                       10218, 10278, 23316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 37941, 0, 3,
                                                                       36621, 22536, 36711,
                                                                       10278, 10338, 23416,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 38091, 0, 3,
                                                                       36711, 22596, 36801,
                                                                       10338, 10398, 23516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 38241, 0, 3,
                                                                       36891, 22716, 37041,
                                                                       10518, 10608, 23616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 38466, 0, 3,
                                                                       37041, 22816, 37191,
                                                                       10608, 10698, 23766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 38691, 0, 3,
                                                                       37191, 22916, 37341,
                                                                       10698, 10788, 23916,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 38916, 0, 3,
                                                                       37341, 23016, 37491,
                                                                       10788, 10878, 24066,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 39141, 0, 3,
                                                                       37491, 23116, 37641,
                                                                       10878, 10968, 24216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 39366, 0, 3,
                                                                       37641, 23216, 37791,
                                                                       10968, 11058, 24366,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 39591, 0, 3,
                                                                       37791, 23316, 37941,
                                                                       11058, 11148, 24516,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 39816, 0, 3,
                                                                       37941, 23416, 38091,
                                                                       11148, 11238, 24666,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 40041, 0, 3,
                                                                       38241, 23616, 38466,
                                                                       11418, 11544, 24816,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 40356, 0, 3,
                                                                       38466, 23766, 38691,
                                                                       11544, 11670, 25026,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 40671, 0, 3,
                                                                       38691, 23916, 38916,
                                                                       11670, 11796, 25236,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 40986, 0, 3,
                                                                       38916, 24066, 39141,
                                                                       11796, 11922, 25446,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 41301, 0, 3,
                                                                       39141, 24216, 39366,
                                                                       11922, 12048, 25656,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 41616, 0, 3,
                                                                       39366, 24366, 39591,
                                                                       12048, 12174, 25866,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 41931, 0, 3,
                                                                       39591, 24516, 39816,
                                                                       12174, 12300, 26076,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 42246, 0, 3,
                                                                       40041, 24816, 40356,
                                                                       12552, 12720, 26286,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 42666, 0, 3,
                                                                       40356, 25026, 40671,
                                                                       12720, 12888, 26566,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 43086, 0, 3,
                                                                       40671, 25236, 40986,
                                                                       12888, 13056, 26846,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 43506, 0, 3,
                                                                       40986, 25446, 41301,
                                                                       13056, 13224, 27126,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 43926, 0, 3,
                                                                       41301, 25656, 41616,
                                                                       13224, 13392, 27406,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 44346, 0, 3,
                                                                       41616, 25866, 41931,
                                                                       13392, 13560, 27686,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 44766, 0, 3,
                                                                       42246, 26286, 42666,
                                                                       13896, 14112, 27966,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 45306, 0, 3,
                                                                       42666, 26566, 43086,
                                                                       14112, 14328, 28326,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 45846, 0, 3,
                                                                       43086, 26846, 43506,
                                                                       14328, 14544, 28686,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 46386, 0, 3,
                                                                       43506, 27126, 43926,
                                                                       14544, 14760, 29046,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 46926, 0, 3,
                                                                       43926, 27406, 44346,
                                                                       14760, 14976, 29406,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 47466, 0, 3,
                                                                       44766, 27966, 45306,
                                                                       15408, 15678, 29766,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 48141, 0, 3,
                                                                       45306, 28326, 45846,
                                                                       15678, 15948, 30216,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 48816, 0, 3,
                                                                       45846, 28686, 46386,
                                                                       15948, 16218, 30666,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 49491, 0, 3,
                                                                       46386, 29046, 46926,
                                                                       16218, 16488, 31116,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 50166, 0, 3,
                                                                       47466, 29766, 48141,
                                                                       17028, 17358, 31566,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 50991, 0, 3,
                                                                       48141, 30216, 48816,
                                                                       17358, 17688, 32116,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 51816, 0, 3,
                                                                       48816, 30666, 49491,
                                                                       17688, 18018, 32666,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 52641, 0, 3,
                                                                       50166, 31566, 50991,
                                                                       18678, 19074, 33216,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 53631, 0, 3,
                                                                       50991, 32116, 51816,
                                                                       19074, 19470, 33876,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 54621, 0, 3,
                                                                       52641, 33216, 53631,
                                                                       20262, 20730, 34536,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 55791, 42246, 420, ncols);

                    simdfunc::contract_primitives(buffer, 56463, 44766, 540, ncols);

                    simdfunc::contract_primitives(buffer, 57327, 47466, 675, ncols);

                    simdfunc::contract_primitives(buffer, 58407, 50166, 825, ncols);

                    simdfunc::contract_primitives(buffer, 59727, 52641, 990, ncols);

                    simdfunc::contract_primitives(buffer, 61311, 54621, 1170, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 56211, 55791, 28, 1, nmax);

        simdtrf::transform_g_inner(buffer, 57003, 56463, 36, 1, nmax);

        simdtrf::transform_g_inner(buffer, 58002, 57327, 45, 1, nmax);

        simdtrf::transform_g_inner(buffer, 59232, 58407, 55, 1, nmax);

        simdtrf::transform_g_inner(buffer, 60717, 59727, 66, 1, nmax);

        simdtrf::transform_g_inner(buffer, 62481, 61311, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 63183, 56211, 57003, 9, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 63939, 57003, 58002, 9, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 64911, 58002, 59232, 9, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 66126, 59232, 60717, 9, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 67611, 60717, 62481, 9, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 69393, 63183, 63939, 9, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 70905, 63939, 64911, 9, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 72849, 64911, 66126, 9, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 75279, 66126, 67611, 9, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 78249, 69393, 70905, 9, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 80769, 70905, 72849, 9, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 84009, 72849, 75279, 9, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 88059, 78249, 80769, 9, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 91839, 80769, 84009, 9, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 96699, 88059, 91839, 9, nmax);

        simdtrf::transform_i_inner(buffer, 101991, 96699, 21, 9, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 101991, 117, nmax);
    }

    for (size_t m = 0; m < 1287; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
