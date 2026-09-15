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


#include "SimdThreeCenterElectronRepulsionRecFIK.hpp"

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
#include "SimdTransferFI.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_fik_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_fik_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 146169, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1365 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 146169, 120990, 7539, dimensions);

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
                                                        16}, ncols, fj, 6, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2043, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2046, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2049, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2052, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2055, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2058, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2061, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2064, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2067, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2070, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2073, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2076, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2079, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2082, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2085, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2088, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2091, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2100, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2109, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2118, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2127, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2136, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2145, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2154, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2163, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2172, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2181, 3, 20, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2190, 3, 21, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2199, 3, 22, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2208, 3, 24, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2226, 3, 27, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2244, 3, 30, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2262, 3, 33, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2280, 3, 36, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2298, 3, 39, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2316, 3, 42, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2334, 3, 45, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2352, 3, 48, 117,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2370, 3, 51, 123,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2388, 3, 54, 129,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2406, 3, 57, 135,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2424, 3, 60, 141,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2442, 3, 63, 147,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2460, 3, 69, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2490, 3, 75, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2520, 3, 81, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2550, 3, 87, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2580, 3, 93, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2610, 3, 99, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2640, 3, 105, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2670, 3, 111, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2700, 3, 117, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2730, 3, 123, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2760, 3, 129, 253,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2790, 3, 135, 263,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2820, 3, 141, 273,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2850, 3, 153, 283,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2895, 3, 163, 298,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2940, 3, 173, 313,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2985, 3, 183, 328,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3030, 3, 193, 343,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3075, 3, 203, 358,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3120, 3, 213, 373,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3165, 3, 223, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3210, 3, 233, 403,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3255, 3, 243, 418,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3300, 3, 253, 433,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3345, 3, 263, 448,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3390, 3, 283, 463,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3453, 3, 298, 484,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3516, 3, 313, 505,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3579, 3, 328, 526,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3642, 3, 343, 547,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3705, 3, 358, 568,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3768, 3, 373, 589,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3831, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3894, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3957, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4020, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4083, 3, 463, 694,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4167, 3, 484, 722,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4251, 3, 505, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4335, 3, 526, 778,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4419, 3, 547, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4503, 3, 568, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4587, 3, 589, 862,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4671, 3, 610, 890,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4755, 3, 631, 918,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4839, 3, 652, 946,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4923, 3, 694, 974,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5031, 3, 722,
                                                                       1010, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5139, 3, 750,
                                                                       1046, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5247, 3, 778,
                                                                       1082, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5355, 3, 806,
                                                                       1118, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5463, 3, 834,
                                                                       1154, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5571, 3, 862,
                                                                       1190, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5679, 3, 890,
                                                                       1226, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5787, 3, 918,
                                                                       1262, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5895, 3, 974,
                                                                       1298, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6030, 3, 1010,
                                                                       1343, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6165, 3, 1046,
                                                                       1388, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6300, 3, 1082,
                                                                       1433, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6435, 3, 1118,
                                                                       1478, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6570, 3, 1154,
                                                                       1523, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6705, 3, 1190,
                                                                       1568, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6840, 3, 1226,
                                                                       1613, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6975, 3, 1298,
                                                                       1658, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7140, 3, 1343,
                                                                       1713, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7305, 3, 1388,
                                                                       1768, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7470, 3, 1433,
                                                                       1823, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7635, 3, 1478,
                                                                       1878, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7800, 3, 1523,
                                                                       1933, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7965, 3, 1568,
                                                                       1988, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8130, 3, 8, 9,
                                                                       2049, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8136, 3, 9, 10,
                                                                       2052, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8142, 3, 10, 11,
                                                                       2055, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8148, 3, 11, 12,
                                                                       2058, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8154, 3, 12, 13,
                                                                       2061, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8160, 3, 13, 14,
                                                                       2064, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8166, 3, 14, 15,
                                                                       2067, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8172, 3, 15, 16,
                                                                       2070, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8178, 3, 16, 17,
                                                                       2073, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8184, 3, 17, 18,
                                                                       2076, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8190, 3, 18, 19,
                                                                       2079, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8196, 3, 19, 20,
                                                                       2082, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8202, 3, 20, 21,
                                                                       2085, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8208, 3, 21, 22,
                                                                       2088, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8214, 0, 3, 8130,
                                                                       2049, 8136, 2091, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8232, 0, 3, 8136,
                                                                       2052, 8142, 2100, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8250, 0, 3, 8142,
                                                                       2055, 8148, 2109, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8268, 0, 3, 8148,
                                                                       2058, 8154, 2118, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8286, 0, 3, 8154,
                                                                       2061, 8160, 2127, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8304, 0, 3, 8160,
                                                                       2064, 8166, 2136, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8322, 0, 3, 8166,
                                                                       2067, 8172, 2145, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8340, 0, 3, 8172,
                                                                       2070, 8178, 2154, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8358, 0, 3, 8178,
                                                                       2073, 8184, 2163, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8376, 0, 3, 8184,
                                                                       2076, 8190, 2172, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8394, 0, 3, 8190,
                                                                       2079, 8196, 2181, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8412, 0, 3, 8196,
                                                                       2082, 8202, 2190, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 8430, 0, 3, 8202,
                                                                       2085, 8208, 2199, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8448, 0, 3, 8214,
                                                                       2091, 8232, 69, 75, 2244,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8484, 0, 3, 8232,
                                                                       2100, 8250, 75, 81, 2262,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8520, 0, 3, 8250,
                                                                       2109, 8268, 81, 87, 2280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8556, 0, 3, 8268,
                                                                       2118, 8286, 87, 93, 2298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8592, 0, 3, 8286,
                                                                       2127, 8304, 93, 99, 2316,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8628, 0, 3, 8304,
                                                                       2136, 8322, 99, 105, 2334,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8664, 0, 3, 8322,
                                                                       2145, 8340, 105, 111,
                                                                       2352, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8700, 0, 3, 8340,
                                                                       2154, 8358, 111, 117,
                                                                       2370, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8736, 0, 3, 8358,
                                                                       2163, 8376, 117, 123,
                                                                       2388, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8772, 0, 3, 8376,
                                                                       2172, 8394, 123, 129,
                                                                       2406, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8808, 0, 3, 8394,
                                                                       2181, 8412, 129, 135,
                                                                       2424, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 8844, 0, 3, 8412,
                                                                       2190, 8430, 135, 141,
                                                                       2442, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8880, 0, 3, 8448,
                                                                       2244, 8484, 153, 163,
                                                                       2520, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8940, 0, 3, 8484,
                                                                       2262, 8520, 163, 173,
                                                                       2550, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9000, 0, 3, 8520,
                                                                       2280, 8556, 173, 183,
                                                                       2580, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9060, 0, 3, 8556,
                                                                       2298, 8592, 183, 193,
                                                                       2610, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9120, 0, 3, 8592,
                                                                       2316, 8628, 193, 203,
                                                                       2640, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9180, 0, 3, 8628,
                                                                       2334, 8664, 203, 213,
                                                                       2670, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9240, 0, 3, 8664,
                                                                       2352, 8700, 213, 223,
                                                                       2700, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9300, 0, 3, 8700,
                                                                       2370, 8736, 223, 233,
                                                                       2730, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9360, 0, 3, 8736,
                                                                       2388, 8772, 233, 243,
                                                                       2760, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9420, 0, 3, 8772,
                                                                       2406, 8808, 243, 253,
                                                                       2790, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 9480, 0, 3, 8808,
                                                                       2424, 8844, 253, 263,
                                                                       2820, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9540, 0, 3, 8880,
                                                                       2520, 8940, 283, 298,
                                                                       2940, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9630, 0, 3, 8940,
                                                                       2550, 9000, 298, 313,
                                                                       2985, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9720, 0, 3, 9000,
                                                                       2580, 9060, 313, 328,
                                                                       3030, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9810, 0, 3, 9060,
                                                                       2610, 9120, 328, 343,
                                                                       3075, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9900, 0, 3, 9120,
                                                                       2640, 9180, 343, 358,
                                                                       3120, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9990, 0, 3, 9180,
                                                                       2670, 9240, 358, 373,
                                                                       3165, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10080, 0, 3, 9240,
                                                                       2700, 9300, 373, 388,
                                                                       3210, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10170, 0, 3, 9300,
                                                                       2730, 9360, 388, 403,
                                                                       3255, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10260, 0, 3, 9360,
                                                                       2760, 9420, 403, 418,
                                                                       3300, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 10350, 0, 3, 9420,
                                                                       2790, 9480, 418, 433,
                                                                       3345, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10440, 0, 3, 9540,
                                                                       2940, 9630, 463, 484,
                                                                       3516, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10566, 0, 3, 9630,
                                                                       2985, 9720, 484, 505,
                                                                       3579, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10692, 0, 3, 9720,
                                                                       3030, 9810, 505, 526,
                                                                       3642, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10818, 0, 3, 9810,
                                                                       3075, 9900, 526, 547,
                                                                       3705, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10944, 0, 3, 9900,
                                                                       3120, 9990, 547, 568,
                                                                       3768, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11070, 0, 3, 9990,
                                                                       3165, 10080, 568, 589,
                                                                       3831, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11196, 0, 3,
                                                                       10080, 3210, 10170, 589,
                                                                       610, 3894, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11322, 0, 3,
                                                                       10170, 3255, 10260, 610,
                                                                       631, 3957, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 11448, 0, 3,
                                                                       10260, 3300, 10350, 631,
                                                                       652, 4020, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11574, 0, 3,
                                                                       10440, 3516, 10566, 694,
                                                                       722, 4251, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11742, 0, 3,
                                                                       10566, 3579, 10692, 722,
                                                                       750, 4335, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11910, 0, 3,
                                                                       10692, 3642, 10818, 750,
                                                                       778, 4419, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12078, 0, 3,
                                                                       10818, 3705, 10944, 778,
                                                                       806, 4503, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12246, 0, 3,
                                                                       10944, 3768, 11070, 806,
                                                                       834, 4587, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12414, 0, 3,
                                                                       11070, 3831, 11196, 834,
                                                                       862, 4671, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12582, 0, 3,
                                                                       11196, 3894, 11322, 862,
                                                                       890, 4755, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 12750, 0, 3,
                                                                       11322, 3957, 11448, 890,
                                                                       918, 4839, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12918, 0, 3,
                                                                       11574, 4251, 11742, 974,
                                                                       1010, 5139, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13134, 0, 3,
                                                                       11742, 4335, 11910, 1010,
                                                                       1046, 5247, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13350, 0, 3,
                                                                       11910, 4419, 12078, 1046,
                                                                       1082, 5355, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13566, 0, 3,
                                                                       12078, 4503, 12246, 1082,
                                                                       1118, 5463, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13782, 0, 3,
                                                                       12246, 4587, 12414, 1118,
                                                                       1154, 5571, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 13998, 0, 3,
                                                                       12414, 4671, 12582, 1154,
                                                                       1190, 5679, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 14214, 0, 3,
                                                                       12582, 4755, 12750, 1190,
                                                                       1226, 5787, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14430, 0, 3,
                                                                       12918, 5139, 13134, 1298,
                                                                       1343, 6165, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14700, 0, 3,
                                                                       13134, 5247, 13350, 1343,
                                                                       1388, 6300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14970, 0, 3,
                                                                       13350, 5355, 13566, 1388,
                                                                       1433, 6435, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15240, 0, 3,
                                                                       13566, 5463, 13782, 1433,
                                                                       1478, 6570, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15510, 0, 3,
                                                                       13782, 5571, 13998, 1478,
                                                                       1523, 6705, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 15780, 0, 3,
                                                                       13998, 5679, 14214, 1523,
                                                                       1568, 6840, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16050, 0, 3,
                                                                       14430, 6165, 14700, 1658,
                                                                       1713, 7305, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16380, 0, 3,
                                                                       14700, 6300, 14970, 1713,
                                                                       1768, 7470, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16710, 0, 3,
                                                                       14970, 6435, 15240, 1768,
                                                                       1823, 7635, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 17040, 0, 3,
                                                                       15240, 6570, 15510, 1823,
                                                                       1878, 7800, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 17370, 0, 3,
                                                                       15510, 6705, 15780, 1878,
                                                                       1933, 7965, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17700, 3, 2043,
                                                                       2046, 8130, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17710, 3, 2046,
                                                                       2049, 8136, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17720, 3, 2049,
                                                                       2052, 8142, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17730, 3, 2052,
                                                                       2055, 8148, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17740, 3, 2055,
                                                                       2058, 8154, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17750, 3, 2058,
                                                                       2061, 8160, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17760, 3, 2061,
                                                                       2064, 8166, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17770, 3, 2064,
                                                                       2067, 8172, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17780, 3, 2067,
                                                                       2070, 8178, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17790, 3, 2070,
                                                                       2073, 8184, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17800, 3, 2073,
                                                                       2076, 8190, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17810, 3, 2076,
                                                                       2079, 8196, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17820, 3, 2079,
                                                                       2082, 8202, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 17830, 3, 2082,
                                                                       2085, 8208, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17840, 0, 3,
                                                                       17700, 8130, 17710, 8214,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17870, 0, 3,
                                                                       17710, 8136, 17720, 8232,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17900, 0, 3,
                                                                       17720, 8142, 17730, 8250,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17930, 0, 3,
                                                                       17730, 8148, 17740, 8268,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17960, 0, 3,
                                                                       17740, 8154, 17750, 8286,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 17990, 0, 3,
                                                                       17750, 8160, 17760, 8304,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18020, 0, 3,
                                                                       17760, 8166, 17770, 8322,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18050, 0, 3,
                                                                       17770, 8172, 17780, 8340,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18080, 0, 3,
                                                                       17780, 8178, 17790, 8358,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18110, 0, 3,
                                                                       17790, 8184, 17800, 8376,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18140, 0, 3,
                                                                       17800, 8190, 17810, 8394,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18170, 0, 3,
                                                                       17810, 8196, 17820, 8412,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18200, 0, 3,
                                                                       17820, 8202, 17830, 8430,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18230, 0, 3,
                                                                       17840, 8214, 17870, 2208,
                                                                       2226, 8448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18290, 0, 3,
                                                                       17870, 8232, 17900, 2226,
                                                                       2244, 8484, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18350, 0, 3,
                                                                       17900, 8250, 17930, 2244,
                                                                       2262, 8520, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18410, 0, 3,
                                                                       17930, 8268, 17960, 2262,
                                                                       2280, 8556, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18470, 0, 3,
                                                                       17960, 8286, 17990, 2280,
                                                                       2298, 8592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18530, 0, 3,
                                                                       17990, 8304, 18020, 2298,
                                                                       2316, 8628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18590, 0, 3,
                                                                       18020, 8322, 18050, 2316,
                                                                       2334, 8664, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18650, 0, 3,
                                                                       18050, 8340, 18080, 2334,
                                                                       2352, 8700, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18710, 0, 3,
                                                                       18080, 8358, 18110, 2352,
                                                                       2370, 8736, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18770, 0, 3,
                                                                       18110, 8376, 18140, 2370,
                                                                       2388, 8772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18830, 0, 3,
                                                                       18140, 8394, 18170, 2388,
                                                                       2406, 8808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 18890, 0, 3,
                                                                       18170, 8412, 18200, 2406,
                                                                       2424, 8844, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18950, 0, 3,
                                                                       18230, 8448, 18290, 2460,
                                                                       2490, 8880, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19050, 0, 3,
                                                                       18290, 8484, 18350, 2490,
                                                                       2520, 8940, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19150, 0, 3,
                                                                       18350, 8520, 18410, 2520,
                                                                       2550, 9000, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19250, 0, 3,
                                                                       18410, 8556, 18470, 2550,
                                                                       2580, 9060, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19350, 0, 3,
                                                                       18470, 8592, 18530, 2580,
                                                                       2610, 9120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19450, 0, 3,
                                                                       18530, 8628, 18590, 2610,
                                                                       2640, 9180, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19550, 0, 3,
                                                                       18590, 8664, 18650, 2640,
                                                                       2670, 9240, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19650, 0, 3,
                                                                       18650, 8700, 18710, 2670,
                                                                       2700, 9300, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19750, 0, 3,
                                                                       18710, 8736, 18770, 2700,
                                                                       2730, 9360, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19850, 0, 3,
                                                                       18770, 8772, 18830, 2730,
                                                                       2760, 9420, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 19950, 0, 3,
                                                                       18830, 8808, 18890, 2760,
                                                                       2790, 9480, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20050, 0, 3,
                                                                       18950, 8880, 19050, 2850,
                                                                       2895, 9540, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20200, 0, 3,
                                                                       19050, 8940, 19150, 2895,
                                                                       2940, 9630, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20350, 0, 3,
                                                                       19150, 9000, 19250, 2940,
                                                                       2985, 9720, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20500, 0, 3,
                                                                       19250, 9060, 19350, 2985,
                                                                       3030, 9810, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20650, 0, 3,
                                                                       19350, 9120, 19450, 3030,
                                                                       3075, 9900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20800, 0, 3,
                                                                       19450, 9180, 19550, 3075,
                                                                       3120, 9990, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 20950, 0, 3,
                                                                       19550, 9240, 19650, 3120,
                                                                       3165, 10080, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21100, 0, 3,
                                                                       19650, 9300, 19750, 3165,
                                                                       3210, 10170, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21250, 0, 3,
                                                                       19750, 9360, 19850, 3210,
                                                                       3255, 10260, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 21400, 0, 3,
                                                                       19850, 9420, 19950, 3255,
                                                                       3300, 10350, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21550, 0, 3,
                                                                       20050, 9540, 20200, 3390,
                                                                       3453, 10440, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21760, 0, 3,
                                                                       20200, 9630, 20350, 3453,
                                                                       3516, 10566, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 21970, 0, 3,
                                                                       20350, 9720, 20500, 3516,
                                                                       3579, 10692, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22180, 0, 3,
                                                                       20500, 9810, 20650, 3579,
                                                                       3642, 10818, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22390, 0, 3,
                                                                       20650, 9900, 20800, 3642,
                                                                       3705, 10944, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22600, 0, 3,
                                                                       20800, 9990, 20950, 3705,
                                                                       3768, 11070, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 22810, 0, 3,
                                                                       20950, 10080, 21100, 3768,
                                                                       3831, 11196, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 23020, 0, 3,
                                                                       21100, 10170, 21250, 3831,
                                                                       3894, 11322, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 23230, 0, 3,
                                                                       21250, 10260, 21400, 3894,
                                                                       3957, 11448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23440, 0, 3,
                                                                       21550, 10440, 21760, 4083,
                                                                       4167, 11574, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 23720, 0, 3,
                                                                       21760, 10566, 21970, 4167,
                                                                       4251, 11742, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24000, 0, 3,
                                                                       21970, 10692, 22180, 4251,
                                                                       4335, 11910, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24280, 0, 3,
                                                                       22180, 10818, 22390, 4335,
                                                                       4419, 12078, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24560, 0, 3,
                                                                       22390, 10944, 22600, 4419,
                                                                       4503, 12246, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 24840, 0, 3,
                                                                       22600, 11070, 22810, 4503,
                                                                       4587, 12414, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 25120, 0, 3,
                                                                       22810, 11196, 23020, 4587,
                                                                       4671, 12582, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 25400, 0, 3,
                                                                       23020, 11322, 23230, 4671,
                                                                       4755, 12750, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 25680, 0, 3,
                                                                       23440, 11574, 23720, 4923,
                                                                       5031, 12918, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26040, 0, 3,
                                                                       23720, 11742, 24000, 5031,
                                                                       5139, 13134, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26400, 0, 3,
                                                                       24000, 11910, 24280, 5139,
                                                                       5247, 13350, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 26760, 0, 3,
                                                                       24280, 12078, 24560, 5247,
                                                                       5355, 13566, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 27120, 0, 3,
                                                                       24560, 12246, 24840, 5355,
                                                                       5463, 13782, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 27480, 0, 3,
                                                                       24840, 12414, 25120, 5463,
                                                                       5571, 13998, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 27840, 0, 3,
                                                                       25120, 12582, 25400, 5571,
                                                                       5679, 14214, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 28200, 0, 3,
                                                                       25680, 12918, 26040, 5895,
                                                                       6030, 14430, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 28650, 0, 3,
                                                                       26040, 13134, 26400, 6030,
                                                                       6165, 14700, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 29100, 0, 3,
                                                                       26400, 13350, 26760, 6165,
                                                                       6300, 14970, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 29550, 0, 3,
                                                                       26760, 13566, 27120, 6300,
                                                                       6435, 15240, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 30000, 0, 3,
                                                                       27120, 13782, 27480, 6435,
                                                                       6570, 15510, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 30450, 0, 3,
                                                                       27480, 13998, 27840, 6570,
                                                                       6705, 15780, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 30900, 0, 3,
                                                                       28200, 14430, 28650, 6975,
                                                                       7140, 16050, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 31450, 0, 3,
                                                                       28650, 14700, 29100, 7140,
                                                                       7305, 16380, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 32000, 0, 3,
                                                                       29100, 14970, 29550, 7305,
                                                                       7470, 16710, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 32550, 0, 3,
                                                                       29550, 15240, 30000, 7470,
                                                                       7635, 17040, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 33100, 0, 3,
                                                                       30000, 15510, 30450, 7635,
                                                                       7800, 17370, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33650, 3, 8130,
                                                                       8136, 17720, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33665, 3, 8136,
                                                                       8142, 17730, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33680, 3, 8142,
                                                                       8148, 17740, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33695, 3, 8148,
                                                                       8154, 17750, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33710, 3, 8154,
                                                                       8160, 17760, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33725, 3, 8160,
                                                                       8166, 17770, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33740, 3, 8166,
                                                                       8172, 17780, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33755, 3, 8172,
                                                                       8178, 17790, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33770, 3, 8178,
                                                                       8184, 17800, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33785, 3, 8184,
                                                                       8190, 17810, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33800, 3, 8190,
                                                                       8196, 17820, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33815, 3, 8196,
                                                                       8202, 17830, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33830, 0, 3,
                                                                       33650, 17720, 33665, 8214,
                                                                       8232, 17900, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33875, 0, 3,
                                                                       33665, 17730, 33680, 8232,
                                                                       8250, 17930, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33920, 0, 3,
                                                                       33680, 17740, 33695, 8250,
                                                                       8268, 17960, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 33965, 0, 3,
                                                                       33695, 17750, 33710, 8268,
                                                                       8286, 17990, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34010, 0, 3,
                                                                       33710, 17760, 33725, 8286,
                                                                       8304, 18020, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34055, 0, 3,
                                                                       33725, 17770, 33740, 8304,
                                                                       8322, 18050, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34100, 0, 3,
                                                                       33740, 17780, 33755, 8322,
                                                                       8340, 18080, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34145, 0, 3,
                                                                       33755, 17790, 33770, 8340,
                                                                       8358, 18110, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34190, 0, 3,
                                                                       33770, 17800, 33785, 8358,
                                                                       8376, 18140, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34235, 0, 3,
                                                                       33785, 17810, 33800, 8376,
                                                                       8394, 18170, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 34280, 0, 3,
                                                                       33800, 17820, 33815, 8394,
                                                                       8412, 18200, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34325, 0, 3,
                                                                       33830, 17900, 33875, 8448,
                                                                       8484, 18350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34415, 0, 3,
                                                                       33875, 17930, 33920, 8484,
                                                                       8520, 18410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34505, 0, 3,
                                                                       33920, 17960, 33965, 8520,
                                                                       8556, 18470, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34595, 0, 3,
                                                                       33965, 17990, 34010, 8556,
                                                                       8592, 18530, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34685, 0, 3,
                                                                       34010, 18020, 34055, 8592,
                                                                       8628, 18590, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34775, 0, 3,
                                                                       34055, 18050, 34100, 8628,
                                                                       8664, 18650, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34865, 0, 3,
                                                                       34100, 18080, 34145, 8664,
                                                                       8700, 18710, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 34955, 0, 3,
                                                                       34145, 18110, 34190, 8700,
                                                                       8736, 18770, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35045, 0, 3,
                                                                       34190, 18140, 34235, 8736,
                                                                       8772, 18830, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 35135, 0, 3,
                                                                       34235, 18170, 34280, 8772,
                                                                       8808, 18890, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35225, 0, 3,
                                                                       34325, 18350, 34415, 8880,
                                                                       8940, 19150, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35375, 0, 3,
                                                                       34415, 18410, 34505, 8940,
                                                                       9000, 19250, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35525, 0, 3,
                                                                       34505, 18470, 34595, 9000,
                                                                       9060, 19350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35675, 0, 3,
                                                                       34595, 18530, 34685, 9060,
                                                                       9120, 19450, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35825, 0, 3,
                                                                       34685, 18590, 34775, 9120,
                                                                       9180, 19550, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 35975, 0, 3,
                                                                       34775, 18650, 34865, 9180,
                                                                       9240, 19650, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36125, 0, 3,
                                                                       34865, 18710, 34955, 9240,
                                                                       9300, 19750, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36275, 0, 3,
                                                                       34955, 18770, 35045, 9300,
                                                                       9360, 19850, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 36425, 0, 3,
                                                                       35045, 18830, 35135, 9360,
                                                                       9420, 19950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 36575, 0, 3,
                                                                       35225, 19150, 35375, 9540,
                                                                       9630, 20350, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 36800, 0, 3,
                                                                       35375, 19250, 35525, 9630,
                                                                       9720, 20500, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37025, 0, 3,
                                                                       35525, 19350, 35675, 9720,
                                                                       9810, 20650, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37250, 0, 3,
                                                                       35675, 19450, 35825, 9810,
                                                                       9900, 20800, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37475, 0, 3,
                                                                       35825, 19550, 35975, 9900,
                                                                       9990, 20950, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37700, 0, 3,
                                                                       35975, 19650, 36125, 9990,
                                                                       10080, 21100, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 37925, 0, 3,
                                                                       36125, 19750, 36275,
                                                                       10080, 10170, 21250,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 38150, 0, 3,
                                                                       36275, 19850, 36425,
                                                                       10170, 10260, 21400,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 38375, 0, 3,
                                                                       36575, 20350, 36800,
                                                                       10440, 10566, 21970,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 38690, 0, 3,
                                                                       36800, 20500, 37025,
                                                                       10566, 10692, 22180,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39005, 0, 3,
                                                                       37025, 20650, 37250,
                                                                       10692, 10818, 22390,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39320, 0, 3,
                                                                       37250, 20800, 37475,
                                                                       10818, 10944, 22600,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39635, 0, 3,
                                                                       37475, 20950, 37700,
                                                                       10944, 11070, 22810,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 39950, 0, 3,
                                                                       37700, 21100, 37925,
                                                                       11070, 11196, 23020,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 40265, 0, 3,
                                                                       37925, 21250, 38150,
                                                                       11196, 11322, 23230,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 40580, 0, 3,
                                                                       38375, 21970, 38690,
                                                                       11574, 11742, 24000,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41000, 0, 3,
                                                                       38690, 22180, 39005,
                                                                       11742, 11910, 24280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41420, 0, 3,
                                                                       39005, 22390, 39320,
                                                                       11910, 12078, 24560,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 41840, 0, 3,
                                                                       39320, 22600, 39635,
                                                                       12078, 12246, 24840,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 42260, 0, 3,
                                                                       39635, 22810, 39950,
                                                                       12246, 12414, 25120,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 42680, 0, 3,
                                                                       39950, 23020, 40265,
                                                                       12414, 12582, 25400,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 43100, 0, 3,
                                                                       40580, 24000, 41000,
                                                                       12918, 13134, 26400,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 43640, 0, 3,
                                                                       41000, 24280, 41420,
                                                                       13134, 13350, 26760,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 44180, 0, 3,
                                                                       41420, 24560, 41840,
                                                                       13350, 13566, 27120,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 44720, 0, 3,
                                                                       41840, 24840, 42260,
                                                                       13566, 13782, 27480,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 45260, 0, 3,
                                                                       42260, 25120, 42680,
                                                                       13782, 13998, 27840,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 45800, 0, 3,
                                                                       43100, 26400, 43640,
                                                                       14430, 14700, 29100,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 46475, 0, 3,
                                                                       43640, 26760, 44180,
                                                                       14700, 14970, 29550,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 47150, 0, 3,
                                                                       44180, 27120, 44720,
                                                                       14970, 15240, 30000,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 47825, 0, 3,
                                                                       44720, 27480, 45260,
                                                                       15240, 15510, 30450,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 48500, 0, 3,
                                                                       45800, 29100, 46475,
                                                                       16050, 16380, 32000,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 49325, 0, 3,
                                                                       46475, 29550, 47150,
                                                                       16380, 16710, 32550,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 50150, 0, 3,
                                                                       47150, 30000, 47825,
                                                                       16710, 17040, 33100,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50975, 3, 17700,
                                                                       17710, 33650, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 50996, 3, 17710,
                                                                       17720, 33665, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51017, 3, 17720,
                                                                       17730, 33680, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51038, 3, 17730,
                                                                       17740, 33695, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51059, 3, 17740,
                                                                       17750, 33710, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51080, 3, 17750,
                                                                       17760, 33725, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51101, 3, 17760,
                                                                       17770, 33740, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51122, 3, 17770,
                                                                       17780, 33755, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51143, 3, 17780,
                                                                       17790, 33770, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51164, 3, 17790,
                                                                       17800, 33785, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51185, 3, 17800,
                                                                       17810, 33800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 51206, 3, 17810,
                                                                       17820, 33815, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51227, 0, 3,
                                                                       50975, 33650, 50996,
                                                                       17840, 17870, 33830,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51290, 0, 3,
                                                                       50996, 33665, 51017,
                                                                       17870, 17900, 33875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51353, 0, 3,
                                                                       51017, 33680, 51038,
                                                                       17900, 17930, 33920,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51416, 0, 3,
                                                                       51038, 33695, 51059,
                                                                       17930, 17960, 33965,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51479, 0, 3,
                                                                       51059, 33710, 51080,
                                                                       17960, 17990, 34010,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51542, 0, 3,
                                                                       51080, 33725, 51101,
                                                                       17990, 18020, 34055,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51605, 0, 3,
                                                                       51101, 33740, 51122,
                                                                       18020, 18050, 34100,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51668, 0, 3,
                                                                       51122, 33755, 51143,
                                                                       18050, 18080, 34145,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51731, 0, 3,
                                                                       51143, 33770, 51164,
                                                                       18080, 18110, 34190,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51794, 0, 3,
                                                                       51164, 33785, 51185,
                                                                       18110, 18140, 34235,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 51857, 0, 3,
                                                                       51185, 33800, 51206,
                                                                       18140, 18170, 34280,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 51920, 0, 3,
                                                                       51227, 33830, 51290,
                                                                       18230, 18290, 34325,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52046, 0, 3,
                                                                       51290, 33875, 51353,
                                                                       18290, 18350, 34415,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52172, 0, 3,
                                                                       51353, 33920, 51416,
                                                                       18350, 18410, 34505,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52298, 0, 3,
                                                                       51416, 33965, 51479,
                                                                       18410, 18470, 34595,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52424, 0, 3,
                                                                       51479, 34010, 51542,
                                                                       18470, 18530, 34685,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52550, 0, 3,
                                                                       51542, 34055, 51605,
                                                                       18530, 18590, 34775,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52676, 0, 3,
                                                                       51605, 34100, 51668,
                                                                       18590, 18650, 34865,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52802, 0, 3,
                                                                       51668, 34145, 51731,
                                                                       18650, 18710, 34955,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 52928, 0, 3,
                                                                       51731, 34190, 51794,
                                                                       18710, 18770, 35045,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 53054, 0, 3,
                                                                       51794, 34235, 51857,
                                                                       18770, 18830, 35135,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 53180, 0, 3,
                                                                       51920, 34325, 52046,
                                                                       18950, 19050, 35225,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 53390, 0, 3,
                                                                       52046, 34415, 52172,
                                                                       19050, 19150, 35375,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 53600, 0, 3,
                                                                       52172, 34505, 52298,
                                                                       19150, 19250, 35525,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 53810, 0, 3,
                                                                       52298, 34595, 52424,
                                                                       19250, 19350, 35675,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54020, 0, 3,
                                                                       52424, 34685, 52550,
                                                                       19350, 19450, 35825,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54230, 0, 3,
                                                                       52550, 34775, 52676,
                                                                       19450, 19550, 35975,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54440, 0, 3,
                                                                       52676, 34865, 52802,
                                                                       19550, 19650, 36125,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54650, 0, 3,
                                                                       52802, 34955, 52928,
                                                                       19650, 19750, 36275,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 54860, 0, 3,
                                                                       52928, 35045, 53054,
                                                                       19750, 19850, 36425,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 55070, 0, 3,
                                                                       53180, 35225, 53390,
                                                                       20050, 20200, 36575,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 55385, 0, 3,
                                                                       53390, 35375, 53600,
                                                                       20200, 20350, 36800,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 55700, 0, 3,
                                                                       53600, 35525, 53810,
                                                                       20350, 20500, 37025,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56015, 0, 3,
                                                                       53810, 35675, 54020,
                                                                       20500, 20650, 37250,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56330, 0, 3,
                                                                       54020, 35825, 54230,
                                                                       20650, 20800, 37475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56645, 0, 3,
                                                                       54230, 35975, 54440,
                                                                       20800, 20950, 37700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 56960, 0, 3,
                                                                       54440, 36125, 54650,
                                                                       20950, 21100, 37925,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 57275, 0, 3,
                                                                       54650, 36275, 54860,
                                                                       21100, 21250, 38150,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 57590, 0, 3,
                                                                       55070, 36575, 55385,
                                                                       21550, 21760, 38375,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 58031, 0, 3,
                                                                       55385, 36800, 55700,
                                                                       21760, 21970, 38690,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 58472, 0, 3,
                                                                       55700, 37025, 56015,
                                                                       21970, 22180, 39005,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 58913, 0, 3,
                                                                       56015, 37250, 56330,
                                                                       22180, 22390, 39320,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 59354, 0, 3,
                                                                       56330, 37475, 56645,
                                                                       22390, 22600, 39635,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 59795, 0, 3,
                                                                       56645, 37700, 56960,
                                                                       22600, 22810, 39950,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 60236, 0, 3,
                                                                       56960, 37925, 57275,
                                                                       22810, 23020, 40265,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 60677, 0, 3,
                                                                       57590, 38375, 58031,
                                                                       23440, 23720, 40580,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 61265, 0, 3,
                                                                       58031, 38690, 58472,
                                                                       23720, 24000, 41000,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 61853, 0, 3,
                                                                       58472, 39005, 58913,
                                                                       24000, 24280, 41420,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 62441, 0, 3,
                                                                       58913, 39320, 59354,
                                                                       24280, 24560, 41840,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 63029, 0, 3,
                                                                       59354, 39635, 59795,
                                                                       24560, 24840, 42260,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 63617, 0, 3,
                                                                       59795, 39950, 60236,
                                                                       24840, 25120, 42680,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 64205, 0, 3,
                                                                       60677, 40580, 61265,
                                                                       25680, 26040, 43100,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 64961, 0, 3,
                                                                       61265, 41000, 61853,
                                                                       26040, 26400, 43640,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 65717, 0, 3,
                                                                       61853, 41420, 62441,
                                                                       26400, 26760, 44180,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 66473, 0, 3,
                                                                       62441, 41840, 63029,
                                                                       26760, 27120, 44720,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 67229, 0, 3,
                                                                       63029, 42260, 63617,
                                                                       27120, 27480, 45260,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 67985, 0, 3,
                                                                       64205, 43100, 64961,
                                                                       28200, 28650, 45800,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 68930, 0, 3,
                                                                       64961, 43640, 65717,
                                                                       28650, 29100, 46475,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 69875, 0, 3,
                                                                       65717, 44180, 66473,
                                                                       29100, 29550, 47150,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 70820, 0, 3,
                                                                       66473, 44720, 67229,
                                                                       29550, 30000, 47825,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 71765, 0, 3,
                                                                       67985, 45800, 68930,
                                                                       30900, 31450, 48500,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 72920, 0, 3,
                                                                       68930, 46475, 69875,
                                                                       31450, 32000, 49325,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 74075, 0, 3,
                                                                       69875, 47150, 70820,
                                                                       32000, 32550, 50150,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75230, 3, 33650,
                                                                       33665, 51017, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75258, 3, 33665,
                                                                       33680, 51038, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75286, 3, 33680,
                                                                       33695, 51059, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75314, 3, 33695,
                                                                       33710, 51080, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75342, 3, 33710,
                                                                       33725, 51101, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75370, 3, 33725,
                                                                       33740, 51122, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75398, 3, 33740,
                                                                       33755, 51143, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75426, 3, 33755,
                                                                       33770, 51164, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75454, 3, 33770,
                                                                       33785, 51185, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 75482, 3, 33785,
                                                                       33800, 51206, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75510, 0, 3,
                                                                       75230, 51017, 75258,
                                                                       33830, 33875, 51353,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75594, 0, 3,
                                                                       75258, 51038, 75286,
                                                                       33875, 33920, 51416,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75678, 0, 3,
                                                                       75286, 51059, 75314,
                                                                       33920, 33965, 51479,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75762, 0, 3,
                                                                       75314, 51080, 75342,
                                                                       33965, 34010, 51542,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75846, 0, 3,
                                                                       75342, 51101, 75370,
                                                                       34010, 34055, 51605,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 75930, 0, 3,
                                                                       75370, 51122, 75398,
                                                                       34055, 34100, 51668,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76014, 0, 3,
                                                                       75398, 51143, 75426,
                                                                       34100, 34145, 51731,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76098, 0, 3,
                                                                       75426, 51164, 75454,
                                                                       34145, 34190, 51794,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 76182, 0, 3,
                                                                       75454, 51185, 75482,
                                                                       34190, 34235, 51857,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76266, 0, 3,
                                                                       75510, 51353, 75594,
                                                                       34325, 34415, 52172,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76434, 0, 3,
                                                                       75594, 51416, 75678,
                                                                       34415, 34505, 52298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76602, 0, 3,
                                                                       75678, 51479, 75762,
                                                                       34505, 34595, 52424,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76770, 0, 3,
                                                                       75762, 51542, 75846,
                                                                       34595, 34685, 52550,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 76938, 0, 3,
                                                                       75846, 51605, 75930,
                                                                       34685, 34775, 52676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77106, 0, 3,
                                                                       75930, 51668, 76014,
                                                                       34775, 34865, 52802,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77274, 0, 3,
                                                                       76014, 51731, 76098,
                                                                       34865, 34955, 52928,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 77442, 0, 3,
                                                                       76098, 51794, 76182,
                                                                       34955, 35045, 53054,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 77610, 0, 3,
                                                                       76266, 52172, 76434,
                                                                       35225, 35375, 53600,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 77890, 0, 3,
                                                                       76434, 52298, 76602,
                                                                       35375, 35525, 53810,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 78170, 0, 3,
                                                                       76602, 52424, 76770,
                                                                       35525, 35675, 54020,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 78450, 0, 3,
                                                                       76770, 52550, 76938,
                                                                       35675, 35825, 54230,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 78730, 0, 3,
                                                                       76938, 52676, 77106,
                                                                       35825, 35975, 54440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 79010, 0, 3,
                                                                       77106, 52802, 77274,
                                                                       35975, 36125, 54650,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 79290, 0, 3,
                                                                       77274, 52928, 77442,
                                                                       36125, 36275, 54860,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 79570, 0, 3,
                                                                       77610, 53600, 77890,
                                                                       36575, 36800, 55700,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 79990, 0, 3,
                                                                       77890, 53810, 78170,
                                                                       36800, 37025, 56015,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 80410, 0, 3,
                                                                       78170, 54020, 78450,
                                                                       37025, 37250, 56330,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 80830, 0, 3,
                                                                       78450, 54230, 78730,
                                                                       37250, 37475, 56645,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 81250, 0, 3,
                                                                       78730, 54440, 79010,
                                                                       37475, 37700, 56960,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 81670, 0, 3,
                                                                       79010, 54650, 79290,
                                                                       37700, 37925, 57275,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 82090, 0, 3,
                                                                       79570, 55700, 79990,
                                                                       38375, 38690, 58472,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 82678, 0, 3,
                                                                       79990, 56015, 80410,
                                                                       38690, 39005, 58913,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 83266, 0, 3,
                                                                       80410, 56330, 80830,
                                                                       39005, 39320, 59354,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 83854, 0, 3,
                                                                       80830, 56645, 81250,
                                                                       39320, 39635, 59795,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 84442, 0, 3,
                                                                       81250, 56960, 81670,
                                                                       39635, 39950, 60236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 85030, 0, 3,
                                                                       82090, 58472, 82678,
                                                                       40580, 41000, 61853,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 85814, 0, 3,
                                                                       82678, 58913, 83266,
                                                                       41000, 41420, 62441,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 86598, 0, 3,
                                                                       83266, 59354, 83854,
                                                                       41420, 41840, 63029,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 87382, 0, 3,
                                                                       83854, 59795, 84442,
                                                                       41840, 42260, 63617,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 88166, 0, 3,
                                                                       85030, 61853, 85814,
                                                                       43100, 43640, 65717,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 89174, 0, 3,
                                                                       85814, 62441, 86598,
                                                                       43640, 44180, 66473,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 90182, 0, 3,
                                                                       86598, 63029, 87382,
                                                                       44180, 44720, 67229,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 91190, 0, 3,
                                                                       88166, 65717, 89174,
                                                                       45800, 46475, 69875,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 92450, 0, 3,
                                                                       89174, 66473, 90182,
                                                                       46475, 47150, 70820,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 93710, 0, 3,
                                                                       91190, 69875, 92450,
                                                                       48500, 49325, 74075,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95250, 3, 50975,
                                                                       50996, 75230, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95286, 3, 50996,
                                                                       51017, 75258, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95322, 3, 51017,
                                                                       51038, 75286, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95358, 3, 51038,
                                                                       51059, 75314, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95394, 3, 51059,
                                                                       51080, 75342, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95430, 3, 51080,
                                                                       51101, 75370, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95466, 3, 51101,
                                                                       51122, 75398, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95502, 3, 51122,
                                                                       51143, 75426, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95538, 3, 51143,
                                                                       51164, 75454, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 95574, 3, 51164,
                                                                       51185, 75482, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95610, 0, 3,
                                                                       95250, 75230, 95286,
                                                                       51227, 51290, 75510,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95718, 0, 3,
                                                                       95286, 75258, 95322,
                                                                       51290, 51353, 75594,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95826, 0, 3,
                                                                       95322, 75286, 95358,
                                                                       51353, 51416, 75678,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 95934, 0, 3,
                                                                       95358, 75314, 95394,
                                                                       51416, 51479, 75762,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 96042, 0, 3,
                                                                       95394, 75342, 95430,
                                                                       51479, 51542, 75846,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 96150, 0, 3,
                                                                       95430, 75370, 95466,
                                                                       51542, 51605, 75930,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 96258, 0, 3,
                                                                       95466, 75398, 95502,
                                                                       51605, 51668, 76014,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 96366, 0, 3,
                                                                       95502, 75426, 95538,
                                                                       51668, 51731, 76098,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 96474, 0, 3,
                                                                       95538, 75454, 95574,
                                                                       51731, 51794, 76182,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 96582, 0, 3,
                                                                       95610, 75510, 95718,
                                                                       51920, 52046, 76266,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 96798, 0, 3,
                                                                       95718, 75594, 95826,
                                                                       52046, 52172, 76434,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97014, 0, 3,
                                                                       95826, 75678, 95934,
                                                                       52172, 52298, 76602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97230, 0, 3,
                                                                       95934, 75762, 96042,
                                                                       52298, 52424, 76770,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97446, 0, 3,
                                                                       96042, 75846, 96150,
                                                                       52424, 52550, 76938,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97662, 0, 3,
                                                                       96150, 75930, 96258,
                                                                       52550, 52676, 77106,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 97878, 0, 3,
                                                                       96258, 76014, 96366,
                                                                       52676, 52802, 77274,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 98094, 0, 3,
                                                                       96366, 76098, 96474,
                                                                       52802, 52928, 77442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 98310, 0, 3,
                                                                       96582, 76266, 96798,
                                                                       53180, 53390, 77610,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 98670, 0, 3,
                                                                       96798, 76434, 97014,
                                                                       53390, 53600, 77890,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 99030, 0, 3,
                                                                       97014, 76602, 97230,
                                                                       53600, 53810, 78170,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 99390, 0, 3,
                                                                       97230, 76770, 97446,
                                                                       53810, 54020, 78450,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 99750, 0, 3,
                                                                       97446, 76938, 97662,
                                                                       54020, 54230, 78730,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 100110, 0, 3,
                                                                       97662, 77106, 97878,
                                                                       54230, 54440, 79010,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 100470, 0, 3,
                                                                       97878, 77274, 98094,
                                                                       54440, 54650, 79290,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 100830, 0, 3,
                                                                       98310, 77610, 98670,
                                                                       55070, 55385, 79570,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 101370, 0, 3,
                                                                       98670, 77890, 99030,
                                                                       55385, 55700, 79990,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 101910, 0, 3,
                                                                       99030, 78170, 99390,
                                                                       55700, 56015, 80410,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 102450, 0, 3,
                                                                       99390, 78450, 99750,
                                                                       56015, 56330, 80830,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 102990, 0, 3,
                                                                       99750, 78730, 100110,
                                                                       56330, 56645, 81250,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 103530, 0, 3,
                                                                       100110, 79010, 100470,
                                                                       56645, 56960, 81670,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 104070, 0, 3,
                                                                       100830, 79570, 101370,
                                                                       57590, 58031, 82090,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 104826, 0, 3,
                                                                       101370, 79990, 101910,
                                                                       58031, 58472, 82678,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 105582, 0, 3,
                                                                       101910, 80410, 102450,
                                                                       58472, 58913, 83266,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 106338, 0, 3,
                                                                       102450, 80830, 102990,
                                                                       58913, 59354, 83854,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 107094, 0, 3,
                                                                       102990, 81250, 103530,
                                                                       59354, 59795, 84442,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 107850, 0, 3,
                                                                       104070, 82090, 104826,
                                                                       60677, 61265, 85030,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 108858, 0, 3,
                                                                       104826, 82678, 105582,
                                                                       61265, 61853, 85814,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 109866, 0, 3,
                                                                       105582, 83266, 106338,
                                                                       61853, 62441, 86598,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 110874, 0, 3,
                                                                       106338, 83854, 107094,
                                                                       62441, 63029, 87382,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 111882, 0, 3,
                                                                       107850, 85030, 108858,
                                                                       64205, 64961, 88166,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 113178, 0, 3,
                                                                       108858, 85814, 109866,
                                                                       64961, 65717, 89174,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 114474, 0, 3,
                                                                       109866, 86598, 110874,
                                                                       65717, 66473, 90182,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 115770, 0, 3,
                                                                       111882, 88166, 113178,
                                                                       67985, 68930, 91190,
                                                                       ncols, gamma, p, q);

                    compute_prim_slk_three_center_electron_repulsion_0(buffer, 117390, 0, 3,
                                                                       113178, 89174, 114474,
                                                                       68930, 69875, 92450,
                                                                       ncols, gamma, p, q);

                    compute_prim_smk_three_center_electron_repulsion_0(buffer, 119010, 0, 3,
                                                                       115770, 91190, 117390,
                                                                       71765, 72920, 93710,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 120990, 107850, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 122418, 111882, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 124254, 115770, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 126549, 119010, 1980, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 121998, 120990, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 123714, 122418, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 125874, 124254, 45, 1, nmax);

        simdtrf::transform_k_inner(buffer, 128529, 126549, 55, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 129354, 121998, 123714, 15, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 130614, 123714, 125874, 15, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 132234, 125874, 128529, 15, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 134259, 129354, 130614, 15, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 136779, 130614, 132234, 15, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 140019, 134259, 136779, 15, nmax);

        simdtrf::transform_i_inner(buffer, 144219, 140019, 10, 15, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 144219, 195, nmax);
    }

    for (size_t m = 0; m < 1365; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
