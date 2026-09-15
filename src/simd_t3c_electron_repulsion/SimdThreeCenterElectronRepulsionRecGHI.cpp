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


#include "SimdThreeCenterElectronRepulsionRecGHI.hpp"

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
#include "SimdTransferDH.hpp"
#include "SimdTransferDI.hpp"
#include "SimdTransferDK.hpp"
#include "SimdTransferFH.hpp"
#include "SimdTransferFI.hpp"
#include "SimdTransferGH.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransferPK.hpp"
#include "SimdTransferPL.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ghi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ghi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 112203, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 112203, 80308, 6870, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2043, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2046, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2049, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2052, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2055, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2058, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2061, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2064, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2067, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2070, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2073, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2076, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2079, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2082, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2085, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2094, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2103, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2112, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2121, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2130, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2139, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2148, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2157, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2166, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2175, 3, 20, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2184, 3, 21, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2193, 3, 22, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2202, 3, 30, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2220, 3, 33, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2238, 3, 36, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2256, 3, 39, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2274, 3, 42, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2292, 3, 45, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2310, 3, 48, 117,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2328, 3, 51, 123,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2346, 3, 54, 129,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2364, 3, 57, 135,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2382, 3, 60, 141,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2400, 3, 63, 147,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2418, 3, 81, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2448, 3, 87, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2478, 3, 93, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2508, 3, 99, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2538, 3, 105, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2568, 3, 111, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2598, 3, 117, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2628, 3, 123, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2658, 3, 129, 253,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2688, 3, 135, 263,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2718, 3, 141, 273,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2748, 3, 173, 313,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2793, 3, 183, 328,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2838, 3, 193, 343,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2883, 3, 203, 358,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2928, 3, 213, 373,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2973, 3, 223, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3018, 3, 233, 403,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3063, 3, 243, 418,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3108, 3, 253, 433,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3153, 3, 263, 448,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3198, 3, 313, 505,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3261, 3, 328, 526,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3324, 3, 343, 547,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3387, 3, 358, 568,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3450, 3, 373, 589,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3513, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3576, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3639, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 3702, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3765, 3, 505, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3849, 3, 526, 778,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3933, 3, 547, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4017, 3, 568, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4101, 3, 589, 862,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4185, 3, 610, 890,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4269, 3, 631, 918,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4353, 3, 652, 946,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4437, 3, 750,
                                                                       1046, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4545, 3, 778,
                                                                       1082, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4653, 3, 806,
                                                                       1118, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4761, 3, 834,
                                                                       1154, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4869, 3, 862,
                                                                       1190, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4977, 3, 890,
                                                                       1226, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5085, 3, 918,
                                                                       1262, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5193, 3, 1046,
                                                                       1388, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5328, 3, 1082,
                                                                       1433, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5463, 3, 1118,
                                                                       1478, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5598, 3, 1154,
                                                                       1523, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5733, 3, 1190,
                                                                       1568, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 5868, 3, 1226,
                                                                       1613, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6003, 3, 1388,
                                                                       1768, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6168, 3, 1433,
                                                                       1823, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6333, 3, 1478,
                                                                       1878, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6498, 3, 1523,
                                                                       1933, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 6663, 3, 1568,
                                                                       1988, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6828, 3, 8, 9,
                                                                       2043, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6834, 3, 9, 10,
                                                                       2046, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6840, 3, 10, 11,
                                                                       2049, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6846, 3, 11, 12,
                                                                       2052, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6852, 3, 12, 13,
                                                                       2055, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6858, 3, 13, 14,
                                                                       2058, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6864, 3, 14, 15,
                                                                       2061, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6870, 3, 15, 16,
                                                                       2064, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6876, 3, 16, 17,
                                                                       2067, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6882, 3, 17, 18,
                                                                       2070, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6888, 3, 18, 19,
                                                                       2073, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6894, 3, 19, 20,
                                                                       2076, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6900, 3, 20, 21,
                                                                       2079, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6906, 3, 21, 22,
                                                                       2082, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6912, 0, 3, 6828,
                                                                       2043, 6834, 2085, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6930, 0, 3, 6834,
                                                                       2046, 6840, 2094, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6948, 0, 3, 6840,
                                                                       2049, 6846, 2103, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6966, 0, 3, 6846,
                                                                       2052, 6852, 2112, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 6984, 0, 3, 6852,
                                                                       2055, 6858, 2121, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7002, 0, 3, 6858,
                                                                       2058, 6864, 2130, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7020, 0, 3, 6864,
                                                                       2061, 6870, 2139, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7038, 0, 3, 6870,
                                                                       2064, 6876, 2148, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7056, 0, 3, 6876,
                                                                       2067, 6882, 2157, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7074, 0, 3, 6882,
                                                                       2070, 6888, 2166, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7092, 0, 3, 6888,
                                                                       2073, 6894, 2175, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7110, 0, 3, 6894,
                                                                       2076, 6900, 2184, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7128, 0, 3, 6900,
                                                                       2079, 6906, 2193, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7146, 0, 3, 6912,
                                                                       2085, 6930, 69, 75, 2202,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7182, 0, 3, 6930,
                                                                       2094, 6948, 75, 81, 2220,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7218, 0, 3, 6948,
                                                                       2103, 6966, 81, 87, 2238,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7254, 0, 3, 6966,
                                                                       2112, 6984, 87, 93, 2256,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7290, 0, 3, 6984,
                                                                       2121, 7002, 93, 99, 2274,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7326, 0, 3, 7002,
                                                                       2130, 7020, 99, 105, 2292,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7362, 0, 3, 7020,
                                                                       2139, 7038, 105, 111,
                                                                       2310, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7398, 0, 3, 7038,
                                                                       2148, 7056, 111, 117,
                                                                       2328, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7434, 0, 3, 7056,
                                                                       2157, 7074, 117, 123,
                                                                       2346, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7470, 0, 3, 7074,
                                                                       2166, 7092, 123, 129,
                                                                       2364, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7506, 0, 3, 7092,
                                                                       2175, 7110, 129, 135,
                                                                       2382, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7542, 0, 3, 7110,
                                                                       2184, 7128, 135, 141,
                                                                       2400, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7578, 0, 3, 7146,
                                                                       2202, 7182, 153, 163,
                                                                       2418, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7638, 0, 3, 7182,
                                                                       2220, 7218, 163, 173,
                                                                       2448, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7698, 0, 3, 7218,
                                                                       2238, 7254, 173, 183,
                                                                       2478, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7758, 0, 3, 7254,
                                                                       2256, 7290, 183, 193,
                                                                       2508, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7818, 0, 3, 7290,
                                                                       2274, 7326, 193, 203,
                                                                       2538, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7878, 0, 3, 7326,
                                                                       2292, 7362, 203, 213,
                                                                       2568, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7938, 0, 3, 7362,
                                                                       2310, 7398, 213, 223,
                                                                       2598, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7998, 0, 3, 7398,
                                                                       2328, 7434, 223, 233,
                                                                       2628, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8058, 0, 3, 7434,
                                                                       2346, 7470, 233, 243,
                                                                       2658, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8118, 0, 3, 7470,
                                                                       2364, 7506, 243, 253,
                                                                       2688, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8178, 0, 3, 7506,
                                                                       2382, 7542, 253, 263,
                                                                       2718, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8238, 0, 3, 7578,
                                                                       2418, 7638, 283, 298,
                                                                       2748, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8328, 0, 3, 7638,
                                                                       2448, 7698, 298, 313,
                                                                       2793, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8418, 0, 3, 7698,
                                                                       2478, 7758, 313, 328,
                                                                       2838, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8508, 0, 3, 7758,
                                                                       2508, 7818, 328, 343,
                                                                       2883, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8598, 0, 3, 7818,
                                                                       2538, 7878, 343, 358,
                                                                       2928, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8688, 0, 3, 7878,
                                                                       2568, 7938, 358, 373,
                                                                       2973, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8778, 0, 3, 7938,
                                                                       2598, 7998, 373, 388,
                                                                       3018, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8868, 0, 3, 7998,
                                                                       2628, 8058, 388, 403,
                                                                       3063, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8958, 0, 3, 8058,
                                                                       2658, 8118, 403, 418,
                                                                       3108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 9048, 0, 3, 8118,
                                                                       2688, 8178, 418, 433,
                                                                       3153, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9138, 0, 3, 8238,
                                                                       2748, 8328, 463, 484,
                                                                       3198, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9264, 0, 3, 8328,
                                                                       2793, 8418, 484, 505,
                                                                       3261, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9390, 0, 3, 8418,
                                                                       2838, 8508, 505, 526,
                                                                       3324, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9516, 0, 3, 8508,
                                                                       2883, 8598, 526, 547,
                                                                       3387, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9642, 0, 3, 8598,
                                                                       2928, 8688, 547, 568,
                                                                       3450, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9768, 0, 3, 8688,
                                                                       2973, 8778, 568, 589,
                                                                       3513, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 9894, 0, 3, 8778,
                                                                       3018, 8868, 589, 610,
                                                                       3576, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10020, 0, 3, 8868,
                                                                       3063, 8958, 610, 631,
                                                                       3639, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 10146, 0, 3, 8958,
                                                                       3108, 9048, 631, 652,
                                                                       3702, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10272, 0, 3, 9138,
                                                                       3198, 9264, 694, 722,
                                                                       3765, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10440, 0, 3, 9264,
                                                                       3261, 9390, 722, 750,
                                                                       3849, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10608, 0, 3, 9390,
                                                                       3324, 9516, 750, 778,
                                                                       3933, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10776, 0, 3, 9516,
                                                                       3387, 9642, 778, 806,
                                                                       4017, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 10944, 0, 3, 9642,
                                                                       3450, 9768, 806, 834,
                                                                       4101, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11112, 0, 3, 9768,
                                                                       3513, 9894, 834, 862,
                                                                       4185, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11280, 0, 3, 9894,
                                                                       3576, 10020, 862, 890,
                                                                       4269, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 11448, 0, 3,
                                                                       10020, 3639, 10146, 890,
                                                                       918, 4353, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11616, 0, 3,
                                                                       10272, 3765, 10440, 974,
                                                                       1010, 4437, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 11832, 0, 3,
                                                                       10440, 3849, 10608, 1010,
                                                                       1046, 4545, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12048, 0, 3,
                                                                       10608, 3933, 10776, 1046,
                                                                       1082, 4653, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12264, 0, 3,
                                                                       10776, 4017, 10944, 1082,
                                                                       1118, 4761, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12480, 0, 3,
                                                                       10944, 4101, 11112, 1118,
                                                                       1154, 4869, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12696, 0, 3,
                                                                       11112, 4185, 11280, 1154,
                                                                       1190, 4977, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 12912, 0, 3,
                                                                       11280, 4269, 11448, 1190,
                                                                       1226, 5085, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13128, 0, 3,
                                                                       11616, 4437, 11832, 1298,
                                                                       1343, 5193, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13398, 0, 3,
                                                                       11832, 4545, 12048, 1343,
                                                                       1388, 5328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13668, 0, 3,
                                                                       12048, 4653, 12264, 1388,
                                                                       1433, 5463, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 13938, 0, 3,
                                                                       12264, 4761, 12480, 1433,
                                                                       1478, 5598, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14208, 0, 3,
                                                                       12480, 4869, 12696, 1478,
                                                                       1523, 5733, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 14478, 0, 3,
                                                                       12696, 4977, 12912, 1523,
                                                                       1568, 5868, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 14748, 0, 3,
                                                                       13128, 5193, 13398, 1658,
                                                                       1713, 6003, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 15078, 0, 3,
                                                                       13398, 5328, 13668, 1713,
                                                                       1768, 6168, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 15408, 0, 3,
                                                                       13668, 5463, 13938, 1768,
                                                                       1823, 6333, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 15738, 0, 3,
                                                                       13938, 5598, 14208, 1823,
                                                                       1878, 6498, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 16068, 0, 3,
                                                                       14208, 5733, 14478, 1878,
                                                                       1933, 6663, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16398, 3, 2043,
                                                                       2046, 6840, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16408, 3, 2046,
                                                                       2049, 6846, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16418, 3, 2049,
                                                                       2052, 6852, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16428, 3, 2052,
                                                                       2055, 6858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16438, 3, 2055,
                                                                       2058, 6864, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16448, 3, 2058,
                                                                       2061, 6870, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16458, 3, 2061,
                                                                       2064, 6876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16468, 3, 2064,
                                                                       2067, 6882, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16478, 3, 2067,
                                                                       2070, 6888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16488, 3, 2070,
                                                                       2073, 6894, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16498, 3, 2073,
                                                                       2076, 6900, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16508, 3, 2076,
                                                                       2079, 6906, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16518, 0, 3,
                                                                       16398, 6840, 16408, 6948,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16548, 0, 3,
                                                                       16408, 6846, 16418, 6966,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16578, 0, 3,
                                                                       16418, 6852, 16428, 6984,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16608, 0, 3,
                                                                       16428, 6858, 16438, 7002,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16638, 0, 3,
                                                                       16438, 6864, 16448, 7020,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16668, 0, 3,
                                                                       16448, 6870, 16458, 7038,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16698, 0, 3,
                                                                       16458, 6876, 16468, 7056,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16728, 0, 3,
                                                                       16468, 6882, 16478, 7074,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16758, 0, 3,
                                                                       16478, 6888, 16488, 7092,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16788, 0, 3,
                                                                       16488, 6894, 16498, 7110,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 16818, 0, 3,
                                                                       16498, 6900, 16508, 7128,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16848, 0, 3,
                                                                       16518, 6948, 16548, 2202,
                                                                       2220, 7218, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16908, 0, 3,
                                                                       16548, 6966, 16578, 2220,
                                                                       2238, 7254, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 16968, 0, 3,
                                                                       16578, 6984, 16608, 2238,
                                                                       2256, 7290, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17028, 0, 3,
                                                                       16608, 7002, 16638, 2256,
                                                                       2274, 7326, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17088, 0, 3,
                                                                       16638, 7020, 16668, 2274,
                                                                       2292, 7362, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17148, 0, 3,
                                                                       16668, 7038, 16698, 2292,
                                                                       2310, 7398, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17208, 0, 3,
                                                                       16698, 7056, 16728, 2310,
                                                                       2328, 7434, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17268, 0, 3,
                                                                       16728, 7074, 16758, 2328,
                                                                       2346, 7470, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17328, 0, 3,
                                                                       16758, 7092, 16788, 2346,
                                                                       2364, 7506, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 17388, 0, 3,
                                                                       16788, 7110, 16818, 2364,
                                                                       2382, 7542, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17448, 0, 3,
                                                                       16848, 7218, 16908, 2418,
                                                                       2448, 7698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17548, 0, 3,
                                                                       16908, 7254, 16968, 2448,
                                                                       2478, 7758, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17648, 0, 3,
                                                                       16968, 7290, 17028, 2478,
                                                                       2508, 7818, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17748, 0, 3,
                                                                       17028, 7326, 17088, 2508,
                                                                       2538, 7878, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17848, 0, 3,
                                                                       17088, 7362, 17148, 2538,
                                                                       2568, 7938, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 17948, 0, 3,
                                                                       17148, 7398, 17208, 2568,
                                                                       2598, 7998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18048, 0, 3,
                                                                       17208, 7434, 17268, 2598,
                                                                       2628, 8058, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18148, 0, 3,
                                                                       17268, 7470, 17328, 2628,
                                                                       2658, 8118, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 18248, 0, 3,
                                                                       17328, 7506, 17388, 2658,
                                                                       2688, 8178, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18348, 0, 3,
                                                                       17448, 7698, 17548, 2748,
                                                                       2793, 8418, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18498, 0, 3,
                                                                       17548, 7758, 17648, 2793,
                                                                       2838, 8508, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18648, 0, 3,
                                                                       17648, 7818, 17748, 2838,
                                                                       2883, 8598, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18798, 0, 3,
                                                                       17748, 7878, 17848, 2883,
                                                                       2928, 8688, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 18948, 0, 3,
                                                                       17848, 7938, 17948, 2928,
                                                                       2973, 8778, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19098, 0, 3,
                                                                       17948, 7998, 18048, 2973,
                                                                       3018, 8868, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19248, 0, 3,
                                                                       18048, 8058, 18148, 3018,
                                                                       3063, 8958, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 19398, 0, 3,
                                                                       18148, 8118, 18248, 3063,
                                                                       3108, 9048, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19548, 0, 3,
                                                                       18348, 8418, 18498, 3198,
                                                                       3261, 9390, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19758, 0, 3,
                                                                       18498, 8508, 18648, 3261,
                                                                       3324, 9516, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 19968, 0, 3,
                                                                       18648, 8598, 18798, 3324,
                                                                       3387, 9642, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20178, 0, 3,
                                                                       18798, 8688, 18948, 3387,
                                                                       3450, 9768, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20388, 0, 3,
                                                                       18948, 8778, 19098, 3450,
                                                                       3513, 9894, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20598, 0, 3,
                                                                       19098, 8868, 19248, 3513,
                                                                       3576, 10020, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 20808, 0, 3,
                                                                       19248, 8958, 19398, 3576,
                                                                       3639, 10146, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21018, 0, 3,
                                                                       19548, 9390, 19758, 3765,
                                                                       3849, 10608, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21298, 0, 3,
                                                                       19758, 9516, 19968, 3849,
                                                                       3933, 10776, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21578, 0, 3,
                                                                       19968, 9642, 20178, 3933,
                                                                       4017, 10944, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 21858, 0, 3,
                                                                       20178, 9768, 20388, 4017,
                                                                       4101, 11112, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22138, 0, 3,
                                                                       20388, 9894, 20598, 4101,
                                                                       4185, 11280, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 22418, 0, 3,
                                                                       20598, 10020, 20808, 4185,
                                                                       4269, 11448, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 22698, 0, 3,
                                                                       21018, 10608, 21298, 4437,
                                                                       4545, 12048, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23058, 0, 3,
                                                                       21298, 10776, 21578, 4545,
                                                                       4653, 12264, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23418, 0, 3,
                                                                       21578, 10944, 21858, 4653,
                                                                       4761, 12480, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 23778, 0, 3,
                                                                       21858, 11112, 22138, 4761,
                                                                       4869, 12696, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 24138, 0, 3,
                                                                       22138, 11280, 22418, 4869,
                                                                       4977, 12912, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 24498, 0, 3,
                                                                       22698, 12048, 23058, 5193,
                                                                       5328, 13668, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 24948, 0, 3,
                                                                       23058, 12264, 23418, 5328,
                                                                       5463, 13938, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 25398, 0, 3,
                                                                       23418, 12480, 23778, 5463,
                                                                       5598, 14208, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 25848, 0, 3,
                                                                       23778, 12696, 24138, 5598,
                                                                       5733, 14478, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 26298, 0, 3,
                                                                       24498, 13668, 24948, 6003,
                                                                       6168, 15408, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 26848, 0, 3,
                                                                       24948, 13938, 25398, 6168,
                                                                       6333, 15738, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 27398, 0, 3,
                                                                       25398, 14208, 25848, 6333,
                                                                       6498, 16068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27948, 3, 6828,
                                                                       6834, 16398, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27963, 3, 6834,
                                                                       6840, 16408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27978, 3, 6840,
                                                                       6846, 16418, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 27993, 3, 6846,
                                                                       6852, 16428, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28008, 3, 6852,
                                                                       6858, 16438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28023, 3, 6858,
                                                                       6864, 16448, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28038, 3, 6864,
                                                                       6870, 16458, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28053, 3, 6870,
                                                                       6876, 16468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28068, 3, 6876,
                                                                       6882, 16478, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28083, 3, 6882,
                                                                       6888, 16488, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28098, 3, 6888,
                                                                       6894, 16498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28113, 3, 6894,
                                                                       6900, 16508, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28128, 0, 3,
                                                                       27948, 16398, 27963, 6912,
                                                                       6930, 16518, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28173, 0, 3,
                                                                       27963, 16408, 27978, 6930,
                                                                       6948, 16548, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28218, 0, 3,
                                                                       27978, 16418, 27993, 6948,
                                                                       6966, 16578, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28263, 0, 3,
                                                                       27993, 16428, 28008, 6966,
                                                                       6984, 16608, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28308, 0, 3,
                                                                       28008, 16438, 28023, 6984,
                                                                       7002, 16638, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28353, 0, 3,
                                                                       28023, 16448, 28038, 7002,
                                                                       7020, 16668, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28398, 0, 3,
                                                                       28038, 16458, 28053, 7020,
                                                                       7038, 16698, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28443, 0, 3,
                                                                       28053, 16468, 28068, 7038,
                                                                       7056, 16728, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28488, 0, 3,
                                                                       28068, 16478, 28083, 7056,
                                                                       7074, 16758, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28533, 0, 3,
                                                                       28083, 16488, 28098, 7074,
                                                                       7092, 16788, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 28578, 0, 3,
                                                                       28098, 16498, 28113, 7092,
                                                                       7110, 16818, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28623, 0, 3,
                                                                       28128, 16518, 28173, 7146,
                                                                       7182, 16848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28713, 0, 3,
                                                                       28173, 16548, 28218, 7182,
                                                                       7218, 16908, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28803, 0, 3,
                                                                       28218, 16578, 28263, 7218,
                                                                       7254, 16968, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28893, 0, 3,
                                                                       28263, 16608, 28308, 7254,
                                                                       7290, 17028, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 28983, 0, 3,
                                                                       28308, 16638, 28353, 7290,
                                                                       7326, 17088, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29073, 0, 3,
                                                                       28353, 16668, 28398, 7326,
                                                                       7362, 17148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29163, 0, 3,
                                                                       28398, 16698, 28443, 7362,
                                                                       7398, 17208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29253, 0, 3,
                                                                       28443, 16728, 28488, 7398,
                                                                       7434, 17268, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29343, 0, 3,
                                                                       28488, 16758, 28533, 7434,
                                                                       7470, 17328, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 29433, 0, 3,
                                                                       28533, 16788, 28578, 7470,
                                                                       7506, 17388, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29523, 0, 3,
                                                                       28623, 16848, 28713, 7578,
                                                                       7638, 17448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29673, 0, 3,
                                                                       28713, 16908, 28803, 7638,
                                                                       7698, 17548, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29823, 0, 3,
                                                                       28803, 16968, 28893, 7698,
                                                                       7758, 17648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 29973, 0, 3,
                                                                       28893, 17028, 28983, 7758,
                                                                       7818, 17748, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30123, 0, 3,
                                                                       28983, 17088, 29073, 7818,
                                                                       7878, 17848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30273, 0, 3,
                                                                       29073, 17148, 29163, 7878,
                                                                       7938, 17948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30423, 0, 3,
                                                                       29163, 17208, 29253, 7938,
                                                                       7998, 18048, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30573, 0, 3,
                                                                       29253, 17268, 29343, 7998,
                                                                       8058, 18148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 30723, 0, 3,
                                                                       29343, 17328, 29433, 8058,
                                                                       8118, 18248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 30873, 0, 3,
                                                                       29523, 17448, 29673, 8238,
                                                                       8328, 18348, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31098, 0, 3,
                                                                       29673, 17548, 29823, 8328,
                                                                       8418, 18498, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31323, 0, 3,
                                                                       29823, 17648, 29973, 8418,
                                                                       8508, 18648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31548, 0, 3,
                                                                       29973, 17748, 30123, 8508,
                                                                       8598, 18798, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31773, 0, 3,
                                                                       30123, 17848, 30273, 8598,
                                                                       8688, 18948, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 31998, 0, 3,
                                                                       30273, 17948, 30423, 8688,
                                                                       8778, 19098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32223, 0, 3,
                                                                       30423, 18048, 30573, 8778,
                                                                       8868, 19248, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 32448, 0, 3,
                                                                       30573, 18148, 30723, 8868,
                                                                       8958, 19398, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 32673, 0, 3,
                                                                       30873, 18348, 31098, 9138,
                                                                       9264, 19548, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 32988, 0, 3,
                                                                       31098, 18498, 31323, 9264,
                                                                       9390, 19758, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 33303, 0, 3,
                                                                       31323, 18648, 31548, 9390,
                                                                       9516, 19968, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 33618, 0, 3,
                                                                       31548, 18798, 31773, 9516,
                                                                       9642, 20178, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 33933, 0, 3,
                                                                       31773, 18948, 31998, 9642,
                                                                       9768, 20388, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 34248, 0, 3,
                                                                       31998, 19098, 32223, 9768,
                                                                       9894, 20598, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 34563, 0, 3,
                                                                       32223, 19248, 32448, 9894,
                                                                       10020, 20808, ncols,
                                                                       gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 34878, 0, 3,
                                                                       32673, 19548, 32988,
                                                                       10272, 10440, 21018,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 35298, 0, 3,
                                                                       32988, 19758, 33303,
                                                                       10440, 10608, 21298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 35718, 0, 3,
                                                                       33303, 19968, 33618,
                                                                       10608, 10776, 21578,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 36138, 0, 3,
                                                                       33618, 20178, 33933,
                                                                       10776, 10944, 21858,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 36558, 0, 3,
                                                                       33933, 20388, 34248,
                                                                       10944, 11112, 22138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 36978, 0, 3,
                                                                       34248, 20598, 34563,
                                                                       11112, 11280, 22418,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 37398, 0, 3,
                                                                       34878, 21018, 35298,
                                                                       11616, 11832, 22698,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 37938, 0, 3,
                                                                       35298, 21298, 35718,
                                                                       11832, 12048, 23058,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 38478, 0, 3,
                                                                       35718, 21578, 36138,
                                                                       12048, 12264, 23418,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 39018, 0, 3,
                                                                       36138, 21858, 36558,
                                                                       12264, 12480, 23778,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 39558, 0, 3,
                                                                       36558, 22138, 36978,
                                                                       12480, 12696, 24138,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 40098, 0, 3,
                                                                       37398, 22698, 37938,
                                                                       13128, 13398, 24498,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 40773, 0, 3,
                                                                       37938, 23058, 38478,
                                                                       13398, 13668, 24948,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 41448, 0, 3,
                                                                       38478, 23418, 39018,
                                                                       13668, 13938, 25398,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 42123, 0, 3,
                                                                       39018, 23778, 39558,
                                                                       13938, 14208, 25848,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 42798, 0, 3,
                                                                       40098, 24498, 40773,
                                                                       14748, 15078, 26298,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 43623, 0, 3,
                                                                       40773, 24948, 41448,
                                                                       15078, 15408, 26848,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 44448, 0, 3,
                                                                       41448, 25398, 42123,
                                                                       15408, 15738, 27398,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45273, 3, 16398,
                                                                       16408, 27978, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45294, 3, 16408,
                                                                       16418, 27993, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45315, 3, 16418,
                                                                       16428, 28008, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45336, 3, 16428,
                                                                       16438, 28023, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45357, 3, 16438,
                                                                       16448, 28038, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45378, 3, 16448,
                                                                       16458, 28053, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45399, 3, 16458,
                                                                       16468, 28068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45420, 3, 16468,
                                                                       16478, 28083, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45441, 3, 16478,
                                                                       16488, 28098, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45462, 3, 16488,
                                                                       16498, 28113, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45483, 0, 3,
                                                                       45273, 27978, 45294,
                                                                       16518, 16548, 28218,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45546, 0, 3,
                                                                       45294, 27993, 45315,
                                                                       16548, 16578, 28263,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45609, 0, 3,
                                                                       45315, 28008, 45336,
                                                                       16578, 16608, 28308,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45672, 0, 3,
                                                                       45336, 28023, 45357,
                                                                       16608, 16638, 28353,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45735, 0, 3,
                                                                       45357, 28038, 45378,
                                                                       16638, 16668, 28398,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45798, 0, 3,
                                                                       45378, 28053, 45399,
                                                                       16668, 16698, 28443,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45861, 0, 3,
                                                                       45399, 28068, 45420,
                                                                       16698, 16728, 28488,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45924, 0, 3,
                                                                       45420, 28083, 45441,
                                                                       16728, 16758, 28533,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 45987, 0, 3,
                                                                       45441, 28098, 45462,
                                                                       16758, 16788, 28578,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46050, 0, 3,
                                                                       45483, 28218, 45546,
                                                                       16848, 16908, 28803,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46176, 0, 3,
                                                                       45546, 28263, 45609,
                                                                       16908, 16968, 28893,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46302, 0, 3,
                                                                       45609, 28308, 45672,
                                                                       16968, 17028, 28983,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46428, 0, 3,
                                                                       45672, 28353, 45735,
                                                                       17028, 17088, 29073,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46554, 0, 3,
                                                                       45735, 28398, 45798,
                                                                       17088, 17148, 29163,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46680, 0, 3,
                                                                       45798, 28443, 45861,
                                                                       17148, 17208, 29253,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46806, 0, 3,
                                                                       45861, 28488, 45924,
                                                                       17208, 17268, 29343,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 46932, 0, 3,
                                                                       45924, 28533, 45987,
                                                                       17268, 17328, 29433,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47058, 0, 3,
                                                                       46050, 28803, 46176,
                                                                       17448, 17548, 29823,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47268, 0, 3,
                                                                       46176, 28893, 46302,
                                                                       17548, 17648, 29973,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47478, 0, 3,
                                                                       46302, 28983, 46428,
                                                                       17648, 17748, 30123,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47688, 0, 3,
                                                                       46428, 29073, 46554,
                                                                       17748, 17848, 30273,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 47898, 0, 3,
                                                                       46554, 29163, 46680,
                                                                       17848, 17948, 30423,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 48108, 0, 3,
                                                                       46680, 29253, 46806,
                                                                       17948, 18048, 30573,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 48318, 0, 3,
                                                                       46806, 29343, 46932,
                                                                       18048, 18148, 30723,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 48528, 0, 3,
                                                                       47058, 29823, 47268,
                                                                       18348, 18498, 31323,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 48843, 0, 3,
                                                                       47268, 29973, 47478,
                                                                       18498, 18648, 31548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 49158, 0, 3,
                                                                       47478, 30123, 47688,
                                                                       18648, 18798, 31773,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 49473, 0, 3,
                                                                       47688, 30273, 47898,
                                                                       18798, 18948, 31998,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 49788, 0, 3,
                                                                       47898, 30423, 48108,
                                                                       18948, 19098, 32223,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 50103, 0, 3,
                                                                       48108, 30573, 48318,
                                                                       19098, 19248, 32448,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 50418, 0, 3,
                                                                       48528, 31323, 48843,
                                                                       19548, 19758, 33303,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 50859, 0, 3,
                                                                       48843, 31548, 49158,
                                                                       19758, 19968, 33618,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 51300, 0, 3,
                                                                       49158, 31773, 49473,
                                                                       19968, 20178, 33933,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 51741, 0, 3,
                                                                       49473, 31998, 49788,
                                                                       20178, 20388, 34248,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 52182, 0, 3,
                                                                       49788, 32223, 50103,
                                                                       20388, 20598, 34563,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 52623, 0, 3,
                                                                       50418, 33303, 50859,
                                                                       21018, 21298, 35718,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 53211, 0, 3,
                                                                       50859, 33618, 51300,
                                                                       21298, 21578, 36138,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 53799, 0, 3,
                                                                       51300, 33933, 51741,
                                                                       21578, 21858, 36558,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 54387, 0, 3,
                                                                       51741, 34248, 52182,
                                                                       21858, 22138, 36978,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 54975, 0, 3,
                                                                       52623, 35718, 53211,
                                                                       22698, 23058, 38478,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 55731, 0, 3,
                                                                       53211, 36138, 53799,
                                                                       23058, 23418, 39018,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 56487, 0, 3,
                                                                       53799, 36558, 54387,
                                                                       23418, 23778, 39558,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 57243, 0, 3,
                                                                       54975, 38478, 55731,
                                                                       24498, 24948, 41448,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 58188, 0, 3,
                                                                       55731, 39018, 56487,
                                                                       24948, 25398, 42123,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 59133, 0, 3,
                                                                       57243, 41448, 58188,
                                                                       26298, 26848, 44448,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60288, 3, 27948,
                                                                       27963, 45273, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60316, 3, 27963,
                                                                       27978, 45294, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60344, 3, 27978,
                                                                       27993, 45315, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60372, 3, 27993,
                                                                       28008, 45336, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60400, 3, 28008,
                                                                       28023, 45357, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60428, 3, 28023,
                                                                       28038, 45378, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60456, 3, 28038,
                                                                       28053, 45399, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60484, 3, 28053,
                                                                       28068, 45420, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60512, 3, 28068,
                                                                       28083, 45441, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 60540, 3, 28083,
                                                                       28098, 45462, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60568, 0, 3,
                                                                       60288, 45273, 60316,
                                                                       28128, 28173, 45483,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60652, 0, 3,
                                                                       60316, 45294, 60344,
                                                                       28173, 28218, 45546,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60736, 0, 3,
                                                                       60344, 45315, 60372,
                                                                       28218, 28263, 45609,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60820, 0, 3,
                                                                       60372, 45336, 60400,
                                                                       28263, 28308, 45672,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60904, 0, 3,
                                                                       60400, 45357, 60428,
                                                                       28308, 28353, 45735,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 60988, 0, 3,
                                                                       60428, 45378, 60456,
                                                                       28353, 28398, 45798,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61072, 0, 3,
                                                                       60456, 45399, 60484,
                                                                       28398, 28443, 45861,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61156, 0, 3,
                                                                       60484, 45420, 60512,
                                                                       28443, 28488, 45924,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61240, 0, 3,
                                                                       60512, 45441, 60540,
                                                                       28488, 28533, 45987,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61324, 0, 3,
                                                                       60568, 45483, 60652,
                                                                       28623, 28713, 46050,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61492, 0, 3,
                                                                       60652, 45546, 60736,
                                                                       28713, 28803, 46176,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61660, 0, 3,
                                                                       60736, 45609, 60820,
                                                                       28803, 28893, 46302,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61828, 0, 3,
                                                                       60820, 45672, 60904,
                                                                       28893, 28983, 46428,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 61996, 0, 3,
                                                                       60904, 45735, 60988,
                                                                       28983, 29073, 46554,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62164, 0, 3,
                                                                       60988, 45798, 61072,
                                                                       29073, 29163, 46680,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62332, 0, 3,
                                                                       61072, 45861, 61156,
                                                                       29163, 29253, 46806,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 62500, 0, 3,
                                                                       61156, 45924, 61240,
                                                                       29253, 29343, 46932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 62668, 0, 3,
                                                                       61324, 46050, 61492,
                                                                       29523, 29673, 47058,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 62948, 0, 3,
                                                                       61492, 46176, 61660,
                                                                       29673, 29823, 47268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63228, 0, 3,
                                                                       61660, 46302, 61828,
                                                                       29823, 29973, 47478,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63508, 0, 3,
                                                                       61828, 46428, 61996,
                                                                       29973, 30123, 47688,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 63788, 0, 3,
                                                                       61996, 46554, 62164,
                                                                       30123, 30273, 47898,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64068, 0, 3,
                                                                       62164, 46680, 62332,
                                                                       30273, 30423, 48108,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 64348, 0, 3,
                                                                       62332, 46806, 62500,
                                                                       30423, 30573, 48318,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 64628, 0, 3,
                                                                       62668, 47058, 62948,
                                                                       30873, 31098, 48528,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 65048, 0, 3,
                                                                       62948, 47268, 63228,
                                                                       31098, 31323, 48843,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 65468, 0, 3,
                                                                       63228, 47478, 63508,
                                                                       31323, 31548, 49158,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 65888, 0, 3,
                                                                       63508, 47688, 63788,
                                                                       31548, 31773, 49473,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 66308, 0, 3,
                                                                       63788, 47898, 64068,
                                                                       31773, 31998, 49788,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 66728, 0, 3,
                                                                       64068, 48108, 64348,
                                                                       31998, 32223, 50103,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 67148, 0, 3,
                                                                       64628, 48528, 65048,
                                                                       32673, 32988, 50418,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 67736, 0, 3,
                                                                       65048, 48843, 65468,
                                                                       32988, 33303, 50859,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 68324, 0, 3,
                                                                       65468, 49158, 65888,
                                                                       33303, 33618, 51300,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 68912, 0, 3,
                                                                       65888, 49473, 66308,
                                                                       33618, 33933, 51741,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 69500, 0, 3,
                                                                       66308, 49788, 66728,
                                                                       33933, 34248, 52182,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 70088, 0, 3,
                                                                       67148, 50418, 67736,
                                                                       34878, 35298, 52623,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 70872, 0, 3,
                                                                       67736, 50859, 68324,
                                                                       35298, 35718, 53211,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 71656, 0, 3,
                                                                       68324, 51300, 68912,
                                                                       35718, 36138, 53799,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 72440, 0, 3,
                                                                       68912, 51741, 69500,
                                                                       36138, 36558, 54387,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 73224, 0, 3,
                                                                       70088, 52623, 70872,
                                                                       37398, 37938, 54975,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 74232, 0, 3,
                                                                       70872, 53211, 71656,
                                                                       37938, 38478, 55731,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 75240, 0, 3,
                                                                       71656, 53799, 72440,
                                                                       38478, 39018, 56487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 76248, 0, 3,
                                                                       73224, 54975, 74232,
                                                                       40098, 40773, 57243,
                                                                       ncols, gamma, p, q);

                    compute_prim_sli_three_center_electron_repulsion_0(buffer, 77508, 0, 3,
                                                                       74232, 55731, 75240,
                                                                       40773, 41448, 58188,
                                                                       ncols, gamma, p, q);

                    compute_prim_smi_three_center_electron_repulsion_0(buffer, 78768, 0, 3,
                                                                       76248, 57243, 77508,
                                                                       42798, 43623, 59133,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 80308, 67148, 588, ncols);

                    simdfunc::contract_primitives(buffer, 81169, 70088, 784, ncols);

                    simdfunc::contract_primitives(buffer, 82317, 73224, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 83793, 76248, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 85638, 78768, 1540, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 80896, 80308, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 81953, 81169, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 83325, 82317, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 85053, 83793, 45, 1, nmax);

        simdtrf::transform_i_inner(buffer, 87178, 85638, 55, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 87893, 80896, 81953, 13, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 88712, 81953, 83325, 13, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 89804, 83325, 85053, 13, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 91208, 85053, 87178, 13, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 92963, 87893, 88712, 13, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 94601, 88712, 89804, 13, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 96785, 89804, 91208, 13, nmax);

        simdtrf::compute_hrr_fh(buffer, coordinates, 99593, 92963, 94601, 13, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 102323, 94601, 96785, 13, nmax);

        simdtrf::compute_hrr_gh(buffer, coordinates, 105963, 99593, 102323, 13, nmax);

        simdtrf::transform_h_inner(buffer, 110058, 105963, 15, 13, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 110058, 143, nmax);
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
