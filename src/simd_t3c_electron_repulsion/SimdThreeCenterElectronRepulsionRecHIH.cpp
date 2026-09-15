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


#include "SimdThreeCenterElectronRepulsionRecHIH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSID.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSND.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSNS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSOS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
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
compute_hih_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hih_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 154193, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1573 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 154193, 93902, 8998, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2829, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2832, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2835, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2838, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2841, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2844, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2847, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2850, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2853, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2856, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2859, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2862, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2865, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2868, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2871, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2874, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2877, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2886, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2895, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2904, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2913, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2922, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2931, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2940, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2949, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2958, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2967, 3, 20, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2976, 3, 21, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2985, 3, 22, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2994, 3, 24, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3012, 3, 27, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3030, 3, 30, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3048, 3, 33, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3066, 3, 36, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3084, 3, 39, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3102, 3, 42, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3120, 3, 45, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3138, 3, 48, 117,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3156, 3, 51, 123,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3174, 3, 54, 129,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3192, 3, 57, 135,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3210, 3, 60, 141,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 3228, 3, 63, 147,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3246, 3, 69, 153,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3276, 3, 75, 163,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3306, 3, 81, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3336, 3, 87, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3366, 3, 93, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3396, 3, 99, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3426, 3, 105, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3456, 3, 111, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3486, 3, 117, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3516, 3, 123, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3546, 3, 129, 253,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3576, 3, 135, 263,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 3606, 3, 141, 273,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3636, 3, 153, 283,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3681, 3, 163, 298,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3726, 3, 173, 313,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3771, 3, 183, 328,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3816, 3, 193, 343,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3861, 3, 203, 358,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3906, 3, 213, 373,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3951, 3, 223, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3996, 3, 233, 403,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4041, 3, 243, 418,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4086, 3, 253, 433,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 4131, 3, 263, 448,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4176, 3, 283, 463,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4239, 3, 298, 484,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4302, 3, 313, 505,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4365, 3, 328, 526,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4428, 3, 343, 547,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4491, 3, 358, 568,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4554, 3, 373, 589,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4617, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4680, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4743, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 4806, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4869, 3, 463, 694,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 4953, 3, 484, 722,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5037, 3, 505, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5121, 3, 526, 778,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5205, 3, 547, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5289, 3, 568, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5373, 3, 589, 862,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5457, 3, 610, 890,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5541, 3, 631, 918,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 5625, 3, 652, 946,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5709, 3, 694, 974,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5817, 3, 722,
                                                                       1010, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 5925, 3, 750,
                                                                       1046, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6033, 3, 778,
                                                                       1082, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6141, 3, 806,
                                                                       1118, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6249, 3, 834,
                                                                       1154, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6357, 3, 862,
                                                                       1190, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6465, 3, 890,
                                                                       1226, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 6573, 3, 918,
                                                                       1262, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6681, 3, 974,
                                                                       1298, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6816, 3, 1010,
                                                                       1343, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 6951, 3, 1046,
                                                                       1388, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7086, 3, 1082,
                                                                       1433, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7221, 3, 1118,
                                                                       1478, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7356, 3, 1154,
                                                                       1523, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7491, 3, 1190,
                                                                       1568, ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 7626, 3, 1226,
                                                                       1613, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7761, 3, 1298,
                                                                       1658, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 7926, 3, 1343,
                                                                       1713, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8091, 3, 1388,
                                                                       1768, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8256, 3, 1433,
                                                                       1823, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8421, 3, 1478,
                                                                       1878, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8586, 3, 1523,
                                                                       1933, ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 8751, 3, 1568,
                                                                       1988, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 8916, 3, 1658,
                                                                       2043, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9114, 3, 1713,
                                                                       2109, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9312, 3, 1768,
                                                                       2175, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9510, 3, 1823,
                                                                       2241, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9708, 3, 1878,
                                                                       2307, ncols, p, q);

                    compute_prim_snp_three_center_electron_repulsion_0(buffer, 9906, 3, 1933,
                                                                       2373, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 10104, 3, 2043,
                                                                       2439, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 10338, 3, 2109,
                                                                       2517, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 10572, 3, 2175,
                                                                       2595, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 10806, 3, 2241,
                                                                       2673, ncols, p, q);

                    compute_prim_sop_three_center_electron_repulsion_0(buffer, 11040, 3, 2307,
                                                                       2751, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11274, 3, 8, 9,
                                                                       2835, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11280, 3, 9, 10,
                                                                       2838, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11286, 3, 10, 11,
                                                                       2841, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11292, 3, 11, 12,
                                                                       2844, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11298, 3, 12, 13,
                                                                       2847, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11304, 3, 13, 14,
                                                                       2850, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11310, 3, 14, 15,
                                                                       2853, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11316, 3, 15, 16,
                                                                       2856, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11322, 3, 16, 17,
                                                                       2859, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11328, 3, 17, 18,
                                                                       2862, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11334, 3, 18, 19,
                                                                       2865, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11340, 3, 19, 20,
                                                                       2868, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11346, 3, 20, 21,
                                                                       2871, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 11352, 3, 21, 22,
                                                                       2874, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11358, 0, 3,
                                                                       11274, 2835, 11280, 2877,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11376, 0, 3,
                                                                       11280, 2838, 11286, 2886,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11394, 0, 3,
                                                                       11286, 2841, 11292, 2895,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11412, 0, 3,
                                                                       11292, 2844, 11298, 2904,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11430, 0, 3,
                                                                       11298, 2847, 11304, 2913,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11448, 0, 3,
                                                                       11304, 2850, 11310, 2922,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11466, 0, 3,
                                                                       11310, 2853, 11316, 2931,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11484, 0, 3,
                                                                       11316, 2856, 11322, 2940,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11502, 0, 3,
                                                                       11322, 2859, 11328, 2949,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11520, 0, 3,
                                                                       11328, 2862, 11334, 2958,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11538, 0, 3,
                                                                       11334, 2865, 11340, 2967,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11556, 0, 3,
                                                                       11340, 2868, 11346, 2976,
                                                                       ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 11574, 0, 3,
                                                                       11346, 2871, 11352, 2985,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11592, 0, 3,
                                                                       11358, 2877, 11376, 69,
                                                                       75, 3030, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11628, 0, 3,
                                                                       11376, 2886, 11394, 75,
                                                                       81, 3048, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11664, 0, 3,
                                                                       11394, 2895, 11412, 81,
                                                                       87, 3066, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11700, 0, 3,
                                                                       11412, 2904, 11430, 87,
                                                                       93, 3084, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11736, 0, 3,
                                                                       11430, 2913, 11448, 93,
                                                                       99, 3102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11772, 0, 3,
                                                                       11448, 2922, 11466, 99,
                                                                       105, 3120, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11808, 0, 3,
                                                                       11466, 2931, 11484, 105,
                                                                       111, 3138, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11844, 0, 3,
                                                                       11484, 2940, 11502, 111,
                                                                       117, 3156, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11880, 0, 3,
                                                                       11502, 2949, 11520, 117,
                                                                       123, 3174, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11916, 0, 3,
                                                                       11520, 2958, 11538, 123,
                                                                       129, 3192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11952, 0, 3,
                                                                       11538, 2967, 11556, 129,
                                                                       135, 3210, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 11988, 0, 3,
                                                                       11556, 2976, 11574, 135,
                                                                       141, 3228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12024, 0, 3,
                                                                       11592, 3030, 11628, 153,
                                                                       163, 3306, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12084, 0, 3,
                                                                       11628, 3048, 11664, 163,
                                                                       173, 3336, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12144, 0, 3,
                                                                       11664, 3066, 11700, 173,
                                                                       183, 3366, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12204, 0, 3,
                                                                       11700, 3084, 11736, 183,
                                                                       193, 3396, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12264, 0, 3,
                                                                       11736, 3102, 11772, 193,
                                                                       203, 3426, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12324, 0, 3,
                                                                       11772, 3120, 11808, 203,
                                                                       213, 3456, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12384, 0, 3,
                                                                       11808, 3138, 11844, 213,
                                                                       223, 3486, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12444, 0, 3,
                                                                       11844, 3156, 11880, 223,
                                                                       233, 3516, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12504, 0, 3,
                                                                       11880, 3174, 11916, 233,
                                                                       243, 3546, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12564, 0, 3,
                                                                       11916, 3192, 11952, 243,
                                                                       253, 3576, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 12624, 0, 3,
                                                                       11952, 3210, 11988, 253,
                                                                       263, 3606, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12684, 0, 3,
                                                                       12024, 3306, 12084, 283,
                                                                       298, 3726, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12774, 0, 3,
                                                                       12084, 3336, 12144, 298,
                                                                       313, 3771, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12864, 0, 3,
                                                                       12144, 3366, 12204, 313,
                                                                       328, 3816, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 12954, 0, 3,
                                                                       12204, 3396, 12264, 328,
                                                                       343, 3861, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13044, 0, 3,
                                                                       12264, 3426, 12324, 343,
                                                                       358, 3906, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13134, 0, 3,
                                                                       12324, 3456, 12384, 358,
                                                                       373, 3951, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13224, 0, 3,
                                                                       12384, 3486, 12444, 373,
                                                                       388, 3996, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13314, 0, 3,
                                                                       12444, 3516, 12504, 388,
                                                                       403, 4041, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13404, 0, 3,
                                                                       12504, 3546, 12564, 403,
                                                                       418, 4086, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 13494, 0, 3,
                                                                       12564, 3576, 12624, 418,
                                                                       433, 4131, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13584, 0, 3,
                                                                       12684, 3726, 12774, 463,
                                                                       484, 4302, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13710, 0, 3,
                                                                       12774, 3771, 12864, 484,
                                                                       505, 4365, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13836, 0, 3,
                                                                       12864, 3816, 12954, 505,
                                                                       526, 4428, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 13962, 0, 3,
                                                                       12954, 3861, 13044, 526,
                                                                       547, 4491, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14088, 0, 3,
                                                                       13044, 3906, 13134, 547,
                                                                       568, 4554, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14214, 0, 3,
                                                                       13134, 3951, 13224, 568,
                                                                       589, 4617, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14340, 0, 3,
                                                                       13224, 3996, 13314, 589,
                                                                       610, 4680, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14466, 0, 3,
                                                                       13314, 4041, 13404, 610,
                                                                       631, 4743, ncols, gamma,
                                                                       p, q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 14592, 0, 3,
                                                                       13404, 4086, 13494, 631,
                                                                       652, 4806, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14718, 0, 3,
                                                                       13584, 4302, 13710, 694,
                                                                       722, 5037, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 14886, 0, 3,
                                                                       13710, 4365, 13836, 722,
                                                                       750, 5121, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15054, 0, 3,
                                                                       13836, 4428, 13962, 750,
                                                                       778, 5205, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15222, 0, 3,
                                                                       13962, 4491, 14088, 778,
                                                                       806, 5289, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15390, 0, 3,
                                                                       14088, 4554, 14214, 806,
                                                                       834, 5373, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15558, 0, 3,
                                                                       14214, 4617, 14340, 834,
                                                                       862, 5457, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15726, 0, 3,
                                                                       14340, 4680, 14466, 862,
                                                                       890, 5541, ncols, gamma,
                                                                       p, q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 15894, 0, 3,
                                                                       14466, 4743, 14592, 890,
                                                                       918, 5625, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16062, 0, 3,
                                                                       14718, 5037, 14886, 974,
                                                                       1010, 5925, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16278, 0, 3,
                                                                       14886, 5121, 15054, 1010,
                                                                       1046, 6033, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16494, 0, 3,
                                                                       15054, 5205, 15222, 1046,
                                                                       1082, 6141, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16710, 0, 3,
                                                                       15222, 5289, 15390, 1082,
                                                                       1118, 6249, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 16926, 0, 3,
                                                                       15390, 5373, 15558, 1118,
                                                                       1154, 6357, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17142, 0, 3,
                                                                       15558, 5457, 15726, 1154,
                                                                       1190, 6465, ncols, gamma,
                                                                       p, q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 17358, 0, 3,
                                                                       15726, 5541, 15894, 1190,
                                                                       1226, 6573, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17574, 0, 3,
                                                                       16062, 5925, 16278, 1298,
                                                                       1343, 6951, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 17844, 0, 3,
                                                                       16278, 6033, 16494, 1343,
                                                                       1388, 7086, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18114, 0, 3,
                                                                       16494, 6141, 16710, 1388,
                                                                       1433, 7221, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18384, 0, 3,
                                                                       16710, 6249, 16926, 1433,
                                                                       1478, 7356, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18654, 0, 3,
                                                                       16926, 6357, 17142, 1478,
                                                                       1523, 7491, ncols, gamma,
                                                                       p, q);

                    compute_prim_sld_three_center_electron_repulsion_0(buffer, 18924, 0, 3,
                                                                       17142, 6465, 17358, 1523,
                                                                       1568, 7626, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19194, 0, 3,
                                                                       17574, 6951, 17844, 1658,
                                                                       1713, 8091, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19524, 0, 3,
                                                                       17844, 7086, 18114, 1713,
                                                                       1768, 8256, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 19854, 0, 3,
                                                                       18114, 7221, 18384, 1768,
                                                                       1823, 8421, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20184, 0, 3,
                                                                       18384, 7356, 18654, 1823,
                                                                       1878, 8586, ncols, gamma,
                                                                       p, q);

                    compute_prim_smd_three_center_electron_repulsion_0(buffer, 20514, 0, 3,
                                                                       18654, 7491, 18924, 1878,
                                                                       1933, 8751, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 20844, 0, 3,
                                                                       19194, 8091, 19524, 2043,
                                                                       2109, 9312, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 21240, 0, 3,
                                                                       19524, 8256, 19854, 2109,
                                                                       2175, 9510, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 21636, 0, 3,
                                                                       19854, 8421, 20184, 2175,
                                                                       2241, 9708, ncols, gamma,
                                                                       p, q);

                    compute_prim_snd_three_center_electron_repulsion_0(buffer, 22032, 0, 3,
                                                                       20184, 8586, 20514, 2241,
                                                                       2307, 9906, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 22428, 0, 3,
                                                                       20844, 9312, 21240, 2439,
                                                                       2517, 10572, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 22896, 0, 3,
                                                                       21240, 9510, 21636, 2517,
                                                                       2595, 10806, ncols, gamma,
                                                                       p, q);

                    compute_prim_sod_three_center_electron_repulsion_0(buffer, 23364, 0, 3,
                                                                       21636, 9708, 22032, 2595,
                                                                       2673, 11040, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23832, 3, 2829,
                                                                       2832, 11274, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23842, 3, 2832,
                                                                       2835, 11280, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23852, 3, 2835,
                                                                       2838, 11286, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23862, 3, 2838,
                                                                       2841, 11292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23872, 3, 2841,
                                                                       2844, 11298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23882, 3, 2844,
                                                                       2847, 11304, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23892, 3, 2847,
                                                                       2850, 11310, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23902, 3, 2850,
                                                                       2853, 11316, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23912, 3, 2853,
                                                                       2856, 11322, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23922, 3, 2856,
                                                                       2859, 11328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23932, 3, 2859,
                                                                       2862, 11334, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23942, 3, 2862,
                                                                       2865, 11340, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23952, 3, 2865,
                                                                       2868, 11346, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 23962, 3, 2868,
                                                                       2871, 11352, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 23972, 0, 3,
                                                                       23832, 11274, 23842,
                                                                       11358, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24002, 0, 3,
                                                                       23842, 11280, 23852,
                                                                       11376, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24032, 0, 3,
                                                                       23852, 11286, 23862,
                                                                       11394, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24062, 0, 3,
                                                                       23862, 11292, 23872,
                                                                       11412, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24092, 0, 3,
                                                                       23872, 11298, 23882,
                                                                       11430, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24122, 0, 3,
                                                                       23882, 11304, 23892,
                                                                       11448, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24152, 0, 3,
                                                                       23892, 11310, 23902,
                                                                       11466, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24182, 0, 3,
                                                                       23902, 11316, 23912,
                                                                       11484, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24212, 0, 3,
                                                                       23912, 11322, 23922,
                                                                       11502, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24242, 0, 3,
                                                                       23922, 11328, 23932,
                                                                       11520, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24272, 0, 3,
                                                                       23932, 11334, 23942,
                                                                       11538, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24302, 0, 3,
                                                                       23942, 11340, 23952,
                                                                       11556, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 24332, 0, 3,
                                                                       23952, 11346, 23962,
                                                                       11574, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24362, 0, 3,
                                                                       23972, 11358, 24002, 2994,
                                                                       3012, 11592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24422, 0, 3,
                                                                       24002, 11376, 24032, 3012,
                                                                       3030, 11628, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24482, 0, 3,
                                                                       24032, 11394, 24062, 3030,
                                                                       3048, 11664, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24542, 0, 3,
                                                                       24062, 11412, 24092, 3048,
                                                                       3066, 11700, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24602, 0, 3,
                                                                       24092, 11430, 24122, 3066,
                                                                       3084, 11736, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24662, 0, 3,
                                                                       24122, 11448, 24152, 3084,
                                                                       3102, 11772, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24722, 0, 3,
                                                                       24152, 11466, 24182, 3102,
                                                                       3120, 11808, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24782, 0, 3,
                                                                       24182, 11484, 24212, 3120,
                                                                       3138, 11844, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24842, 0, 3,
                                                                       24212, 11502, 24242, 3138,
                                                                       3156, 11880, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24902, 0, 3,
                                                                       24242, 11520, 24272, 3156,
                                                                       3174, 11916, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 24962, 0, 3,
                                                                       24272, 11538, 24302, 3174,
                                                                       3192, 11952, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 25022, 0, 3,
                                                                       24302, 11556, 24332, 3192,
                                                                       3210, 11988, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25082, 0, 3,
                                                                       24362, 11592, 24422, 3246,
                                                                       3276, 12024, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25182, 0, 3,
                                                                       24422, 11628, 24482, 3276,
                                                                       3306, 12084, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25282, 0, 3,
                                                                       24482, 11664, 24542, 3306,
                                                                       3336, 12144, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25382, 0, 3,
                                                                       24542, 11700, 24602, 3336,
                                                                       3366, 12204, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25482, 0, 3,
                                                                       24602, 11736, 24662, 3366,
                                                                       3396, 12264, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25582, 0, 3,
                                                                       24662, 11772, 24722, 3396,
                                                                       3426, 12324, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25682, 0, 3,
                                                                       24722, 11808, 24782, 3426,
                                                                       3456, 12384, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25782, 0, 3,
                                                                       24782, 11844, 24842, 3456,
                                                                       3486, 12444, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25882, 0, 3,
                                                                       24842, 11880, 24902, 3486,
                                                                       3516, 12504, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 25982, 0, 3,
                                                                       24902, 11916, 24962, 3516,
                                                                       3546, 12564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 26082, 0, 3,
                                                                       24962, 11952, 25022, 3546,
                                                                       3576, 12624, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26182, 0, 3,
                                                                       25082, 12024, 25182, 3636,
                                                                       3681, 12684, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26332, 0, 3,
                                                                       25182, 12084, 25282, 3681,
                                                                       3726, 12774, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26482, 0, 3,
                                                                       25282, 12144, 25382, 3726,
                                                                       3771, 12864, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26632, 0, 3,
                                                                       25382, 12204, 25482, 3771,
                                                                       3816, 12954, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26782, 0, 3,
                                                                       25482, 12264, 25582, 3816,
                                                                       3861, 13044, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 26932, 0, 3,
                                                                       25582, 12324, 25682, 3861,
                                                                       3906, 13134, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27082, 0, 3,
                                                                       25682, 12384, 25782, 3906,
                                                                       3951, 13224, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27232, 0, 3,
                                                                       25782, 12444, 25882, 3951,
                                                                       3996, 13314, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27382, 0, 3,
                                                                       25882, 12504, 25982, 3996,
                                                                       4041, 13404, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 27532, 0, 3,
                                                                       25982, 12564, 26082, 4041,
                                                                       4086, 13494, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27682, 0, 3,
                                                                       26182, 12684, 26332, 4176,
                                                                       4239, 13584, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 27892, 0, 3,
                                                                       26332, 12774, 26482, 4239,
                                                                       4302, 13710, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28102, 0, 3,
                                                                       26482, 12864, 26632, 4302,
                                                                       4365, 13836, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28312, 0, 3,
                                                                       26632, 12954, 26782, 4365,
                                                                       4428, 13962, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28522, 0, 3,
                                                                       26782, 13044, 26932, 4428,
                                                                       4491, 14088, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28732, 0, 3,
                                                                       26932, 13134, 27082, 4491,
                                                                       4554, 14214, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 28942, 0, 3,
                                                                       27082, 13224, 27232, 4554,
                                                                       4617, 14340, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29152, 0, 3,
                                                                       27232, 13314, 27382, 4617,
                                                                       4680, 14466, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 29362, 0, 3,
                                                                       27382, 13404, 27532, 4680,
                                                                       4743, 14592, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29572, 0, 3,
                                                                       27682, 13584, 27892, 4869,
                                                                       4953, 14718, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 29852, 0, 3,
                                                                       27892, 13710, 28102, 4953,
                                                                       5037, 14886, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30132, 0, 3,
                                                                       28102, 13836, 28312, 5037,
                                                                       5121, 15054, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30412, 0, 3,
                                                                       28312, 13962, 28522, 5121,
                                                                       5205, 15222, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30692, 0, 3,
                                                                       28522, 14088, 28732, 5205,
                                                                       5289, 15390, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 30972, 0, 3,
                                                                       28732, 14214, 28942, 5289,
                                                                       5373, 15558, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31252, 0, 3,
                                                                       28942, 14340, 29152, 5373,
                                                                       5457, 15726, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 31532, 0, 3,
                                                                       29152, 14466, 29362, 5457,
                                                                       5541, 15894, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 31812, 0, 3,
                                                                       29572, 14718, 29852, 5709,
                                                                       5817, 16062, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32172, 0, 3,
                                                                       29852, 14886, 30132, 5817,
                                                                       5925, 16278, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32532, 0, 3,
                                                                       30132, 15054, 30412, 5925,
                                                                       6033, 16494, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 32892, 0, 3,
                                                                       30412, 15222, 30692, 6033,
                                                                       6141, 16710, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33252, 0, 3,
                                                                       30692, 15390, 30972, 6141,
                                                                       6249, 16926, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33612, 0, 3,
                                                                       30972, 15558, 31252, 6249,
                                                                       6357, 17142, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 33972, 0, 3,
                                                                       31252, 15726, 31532, 6357,
                                                                       6465, 17358, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 34332, 0, 3,
                                                                       31812, 16062, 32172, 6681,
                                                                       6816, 17574, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 34782, 0, 3,
                                                                       32172, 16278, 32532, 6816,
                                                                       6951, 17844, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35232, 0, 3,
                                                                       32532, 16494, 32892, 6951,
                                                                       7086, 18114, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 35682, 0, 3,
                                                                       32892, 16710, 33252, 7086,
                                                                       7221, 18384, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 36132, 0, 3,
                                                                       33252, 16926, 33612, 7221,
                                                                       7356, 18654, ncols, gamma,
                                                                       p, q);

                    compute_prim_slf_three_center_electron_repulsion_0(buffer, 36582, 0, 3,
                                                                       33612, 17142, 33972, 7356,
                                                                       7491, 18924, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 37032, 0, 3,
                                                                       34332, 17574, 34782, 7761,
                                                                       7926, 19194, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 37582, 0, 3,
                                                                       34782, 17844, 35232, 7926,
                                                                       8091, 19524, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 38132, 0, 3,
                                                                       35232, 18114, 35682, 8091,
                                                                       8256, 19854, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 38682, 0, 3,
                                                                       35682, 18384, 36132, 8256,
                                                                       8421, 20184, ncols, gamma,
                                                                       p, q);

                    compute_prim_smf_three_center_electron_repulsion_0(buffer, 39232, 0, 3,
                                                                       36132, 18654, 36582, 8421,
                                                                       8586, 20514, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 39782, 0, 3,
                                                                       37032, 19194, 37582, 8916,
                                                                       9114, 20844, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 40442, 0, 3,
                                                                       37582, 19524, 38132, 9114,
                                                                       9312, 21240, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 41102, 0, 3,
                                                                       38132, 19854, 38682, 9312,
                                                                       9510, 21636, ncols, gamma,
                                                                       p, q);

                    compute_prim_snf_three_center_electron_repulsion_0(buffer, 41762, 0, 3,
                                                                       38682, 20184, 39232, 9510,
                                                                       9708, 22032, ncols, gamma,
                                                                       p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 42422, 0, 3,
                                                                       39782, 20844, 40442,
                                                                       10104, 10338, 22428,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 43202, 0, 3,
                                                                       40442, 21240, 41102,
                                                                       10338, 10572, 22896,
                                                                       ncols, gamma, p, q);

                    compute_prim_sof_three_center_electron_repulsion_0(buffer, 43982, 0, 3,
                                                                       41102, 21636, 41762,
                                                                       10572, 10806, 23364,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44762, 3, 11274,
                                                                       11280, 23852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44777, 3, 11280,
                                                                       11286, 23862, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44792, 3, 11286,
                                                                       11292, 23872, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44807, 3, 11292,
                                                                       11298, 23882, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44822, 3, 11298,
                                                                       11304, 23892, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44837, 3, 11304,
                                                                       11310, 23902, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44852, 3, 11310,
                                                                       11316, 23912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44867, 3, 11316,
                                                                       11322, 23922, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44882, 3, 11322,
                                                                       11328, 23932, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44897, 3, 11328,
                                                                       11334, 23942, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44912, 3, 11334,
                                                                       11340, 23952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 44927, 3, 11340,
                                                                       11346, 23962, ncols,
                                                                       gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 44942, 0, 3,
                                                                       44762, 23852, 44777,
                                                                       11358, 11376, 24032,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 44987, 0, 3,
                                                                       44777, 23862, 44792,
                                                                       11376, 11394, 24062,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45032, 0, 3,
                                                                       44792, 23872, 44807,
                                                                       11394, 11412, 24092,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45077, 0, 3,
                                                                       44807, 23882, 44822,
                                                                       11412, 11430, 24122,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45122, 0, 3,
                                                                       44822, 23892, 44837,
                                                                       11430, 11448, 24152,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45167, 0, 3,
                                                                       44837, 23902, 44852,
                                                                       11448, 11466, 24182,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45212, 0, 3,
                                                                       44852, 23912, 44867,
                                                                       11466, 11484, 24212,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45257, 0, 3,
                                                                       44867, 23922, 44882,
                                                                       11484, 11502, 24242,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45302, 0, 3,
                                                                       44882, 23932, 44897,
                                                                       11502, 11520, 24272,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45347, 0, 3,
                                                                       44897, 23942, 44912,
                                                                       11520, 11538, 24302,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 45392, 0, 3,
                                                                       44912, 23952, 44927,
                                                                       11538, 11556, 24332,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45437, 0, 3,
                                                                       44942, 24032, 44987,
                                                                       11592, 11628, 24482,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45527, 0, 3,
                                                                       44987, 24062, 45032,
                                                                       11628, 11664, 24542,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45617, 0, 3,
                                                                       45032, 24092, 45077,
                                                                       11664, 11700, 24602,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45707, 0, 3,
                                                                       45077, 24122, 45122,
                                                                       11700, 11736, 24662,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45797, 0, 3,
                                                                       45122, 24152, 45167,
                                                                       11736, 11772, 24722,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45887, 0, 3,
                                                                       45167, 24182, 45212,
                                                                       11772, 11808, 24782,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 45977, 0, 3,
                                                                       45212, 24212, 45257,
                                                                       11808, 11844, 24842,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46067, 0, 3,
                                                                       45257, 24242, 45302,
                                                                       11844, 11880, 24902,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46157, 0, 3,
                                                                       45302, 24272, 45347,
                                                                       11880, 11916, 24962,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 46247, 0, 3,
                                                                       45347, 24302, 45392,
                                                                       11916, 11952, 25022,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46337, 0, 3,
                                                                       45437, 24482, 45527,
                                                                       12024, 12084, 25282,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46487, 0, 3,
                                                                       45527, 24542, 45617,
                                                                       12084, 12144, 25382,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46637, 0, 3,
                                                                       45617, 24602, 45707,
                                                                       12144, 12204, 25482,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46787, 0, 3,
                                                                       45707, 24662, 45797,
                                                                       12204, 12264, 25582,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 46937, 0, 3,
                                                                       45797, 24722, 45887,
                                                                       12264, 12324, 25682,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47087, 0, 3,
                                                                       45887, 24782, 45977,
                                                                       12324, 12384, 25782,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47237, 0, 3,
                                                                       45977, 24842, 46067,
                                                                       12384, 12444, 25882,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47387, 0, 3,
                                                                       46067, 24902, 46157,
                                                                       12444, 12504, 25982,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 47537, 0, 3,
                                                                       46157, 24962, 46247,
                                                                       12504, 12564, 26082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 47687, 0, 3,
                                                                       46337, 25282, 46487,
                                                                       12684, 12774, 26482,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 47912, 0, 3,
                                                                       46487, 25382, 46637,
                                                                       12774, 12864, 26632,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48137, 0, 3,
                                                                       46637, 25482, 46787,
                                                                       12864, 12954, 26782,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48362, 0, 3,
                                                                       46787, 25582, 46937,
                                                                       12954, 13044, 26932,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48587, 0, 3,
                                                                       46937, 25682, 47087,
                                                                       13044, 13134, 27082,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 48812, 0, 3,
                                                                       47087, 25782, 47237,
                                                                       13134, 13224, 27232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49037, 0, 3,
                                                                       47237, 25882, 47387,
                                                                       13224, 13314, 27382,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 49262, 0, 3,
                                                                       47387, 25982, 47537,
                                                                       13314, 13404, 27532,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 49487, 0, 3,
                                                                       47687, 26482, 47912,
                                                                       13584, 13710, 28102,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 49802, 0, 3,
                                                                       47912, 26632, 48137,
                                                                       13710, 13836, 28312,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50117, 0, 3,
                                                                       48137, 26782, 48362,
                                                                       13836, 13962, 28522,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50432, 0, 3,
                                                                       48362, 26932, 48587,
                                                                       13962, 14088, 28732,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 50747, 0, 3,
                                                                       48587, 27082, 48812,
                                                                       14088, 14214, 28942,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51062, 0, 3,
                                                                       48812, 27232, 49037,
                                                                       14214, 14340, 29152,
                                                                       ncols, gamma, p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 51377, 0, 3,
                                                                       49037, 27382, 49262,
                                                                       14340, 14466, 29362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 51692, 0, 3,
                                                                       49487, 28102, 49802,
                                                                       14718, 14886, 30132,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 52112, 0, 3,
                                                                       49802, 28312, 50117,
                                                                       14886, 15054, 30412,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 52532, 0, 3,
                                                                       50117, 28522, 50432,
                                                                       15054, 15222, 30692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 52952, 0, 3,
                                                                       50432, 28732, 50747,
                                                                       15222, 15390, 30972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 53372, 0, 3,
                                                                       50747, 28942, 51062,
                                                                       15390, 15558, 31252,
                                                                       ncols, gamma, p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 53792, 0, 3,
                                                                       51062, 29152, 51377,
                                                                       15558, 15726, 31532,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 54212, 0, 3,
                                                                       51692, 30132, 52112,
                                                                       16062, 16278, 32532,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 54752, 0, 3,
                                                                       52112, 30412, 52532,
                                                                       16278, 16494, 32892,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 55292, 0, 3,
                                                                       52532, 30692, 52952,
                                                                       16494, 16710, 33252,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 55832, 0, 3,
                                                                       52952, 30972, 53372,
                                                                       16710, 16926, 33612,
                                                                       ncols, gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 56372, 0, 3,
                                                                       53372, 31252, 53792,
                                                                       16926, 17142, 33972,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 56912, 0, 3,
                                                                       54212, 32532, 54752,
                                                                       17574, 17844, 35232,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 57587, 0, 3,
                                                                       54752, 32892, 55292,
                                                                       17844, 18114, 35682,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 58262, 0, 3,
                                                                       55292, 33252, 55832,
                                                                       18114, 18384, 36132,
                                                                       ncols, gamma, p, q);

                    compute_prim_slg_three_center_electron_repulsion_0(buffer, 58937, 0, 3,
                                                                       55832, 33612, 56372,
                                                                       18384, 18654, 36582,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 59612, 0, 3,
                                                                       56912, 35232, 57587,
                                                                       19194, 19524, 38132,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 60437, 0, 3,
                                                                       57587, 35682, 58262,
                                                                       19524, 19854, 38682,
                                                                       ncols, gamma, p, q);

                    compute_prim_smg_three_center_electron_repulsion_0(buffer, 61262, 0, 3,
                                                                       58262, 36132, 58937,
                                                                       19854, 20184, 39232,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 62087, 0, 3,
                                                                       59612, 38132, 60437,
                                                                       20844, 21240, 41102,
                                                                       ncols, gamma, p, q);

                    compute_prim_sng_three_center_electron_repulsion_0(buffer, 63077, 0, 3,
                                                                       60437, 38682, 61262,
                                                                       21240, 21636, 41762,
                                                                       ncols, gamma, p, q);

                    compute_prim_sog_three_center_electron_repulsion_0(buffer, 64067, 0, 3,
                                                                       62087, 41102, 63077,
                                                                       22428, 22896, 43982,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65237, 3, 23832,
                                                                       23842, 44762, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65258, 3, 23842,
                                                                       23852, 44777, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65279, 3, 23852,
                                                                       23862, 44792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65300, 3, 23862,
                                                                       23872, 44807, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65321, 3, 23872,
                                                                       23882, 44822, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65342, 3, 23882,
                                                                       23892, 44837, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65363, 3, 23892,
                                                                       23902, 44852, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65384, 3, 23902,
                                                                       23912, 44867, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65405, 3, 23912,
                                                                       23922, 44882, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65426, 3, 23922,
                                                                       23932, 44897, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65447, 3, 23932,
                                                                       23942, 44912, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 65468, 3, 23942,
                                                                       23952, 44927, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65489, 0, 3,
                                                                       65237, 44762, 65258,
                                                                       23972, 24002, 44942,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65552, 0, 3,
                                                                       65258, 44777, 65279,
                                                                       24002, 24032, 44987,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65615, 0, 3,
                                                                       65279, 44792, 65300,
                                                                       24032, 24062, 45032,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65678, 0, 3,
                                                                       65300, 44807, 65321,
                                                                       24062, 24092, 45077,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65741, 0, 3,
                                                                       65321, 44822, 65342,
                                                                       24092, 24122, 45122,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65804, 0, 3,
                                                                       65342, 44837, 65363,
                                                                       24122, 24152, 45167,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65867, 0, 3,
                                                                       65363, 44852, 65384,
                                                                       24152, 24182, 45212,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65930, 0, 3,
                                                                       65384, 44867, 65405,
                                                                       24182, 24212, 45257,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 65993, 0, 3,
                                                                       65405, 44882, 65426,
                                                                       24212, 24242, 45302,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 66056, 0, 3,
                                                                       65426, 44897, 65447,
                                                                       24242, 24272, 45347,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 66119, 0, 3,
                                                                       65447, 44912, 65468,
                                                                       24272, 24302, 45392,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66182, 0, 3,
                                                                       65489, 44942, 65552,
                                                                       24362, 24422, 45437,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66308, 0, 3,
                                                                       65552, 44987, 65615,
                                                                       24422, 24482, 45527,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66434, 0, 3,
                                                                       65615, 45032, 65678,
                                                                       24482, 24542, 45617,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66560, 0, 3,
                                                                       65678, 45077, 65741,
                                                                       24542, 24602, 45707,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66686, 0, 3,
                                                                       65741, 45122, 65804,
                                                                       24602, 24662, 45797,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66812, 0, 3,
                                                                       65804, 45167, 65867,
                                                                       24662, 24722, 45887,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 66938, 0, 3,
                                                                       65867, 45212, 65930,
                                                                       24722, 24782, 45977,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 67064, 0, 3,
                                                                       65930, 45257, 65993,
                                                                       24782, 24842, 46067,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 67190, 0, 3,
                                                                       65993, 45302, 66056,
                                                                       24842, 24902, 46157,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 67316, 0, 3,
                                                                       66056, 45347, 66119,
                                                                       24902, 24962, 46247,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67442, 0, 3,
                                                                       66182, 45437, 66308,
                                                                       25082, 25182, 46337,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67652, 0, 3,
                                                                       66308, 45527, 66434,
                                                                       25182, 25282, 46487,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 67862, 0, 3,
                                                                       66434, 45617, 66560,
                                                                       25282, 25382, 46637,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68072, 0, 3,
                                                                       66560, 45707, 66686,
                                                                       25382, 25482, 46787,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68282, 0, 3,
                                                                       66686, 45797, 66812,
                                                                       25482, 25582, 46937,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68492, 0, 3,
                                                                       66812, 45887, 66938,
                                                                       25582, 25682, 47087,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68702, 0, 3,
                                                                       66938, 45977, 67064,
                                                                       25682, 25782, 47237,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 68912, 0, 3,
                                                                       67064, 46067, 67190,
                                                                       25782, 25882, 47387,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 69122, 0, 3,
                                                                       67190, 46157, 67316,
                                                                       25882, 25982, 47537,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 69332, 0, 3,
                                                                       67442, 46337, 67652,
                                                                       26182, 26332, 47687,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 69647, 0, 3,
                                                                       67652, 46487, 67862,
                                                                       26332, 26482, 47912,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 69962, 0, 3,
                                                                       67862, 46637, 68072,
                                                                       26482, 26632, 48137,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 70277, 0, 3,
                                                                       68072, 46787, 68282,
                                                                       26632, 26782, 48362,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 70592, 0, 3,
                                                                       68282, 46937, 68492,
                                                                       26782, 26932, 48587,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 70907, 0, 3,
                                                                       68492, 47087, 68702,
                                                                       26932, 27082, 48812,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 71222, 0, 3,
                                                                       68702, 47237, 68912,
                                                                       27082, 27232, 49037,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 71537, 0, 3,
                                                                       68912, 47387, 69122,
                                                                       27232, 27382, 49262,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 71852, 0, 3,
                                                                       69332, 47687, 69647,
                                                                       27682, 27892, 49487,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 72293, 0, 3,
                                                                       69647, 47912, 69962,
                                                                       27892, 28102, 49802,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 72734, 0, 3,
                                                                       69962, 48137, 70277,
                                                                       28102, 28312, 50117,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 73175, 0, 3,
                                                                       70277, 48362, 70592,
                                                                       28312, 28522, 50432,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 73616, 0, 3,
                                                                       70592, 48587, 70907,
                                                                       28522, 28732, 50747,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 74057, 0, 3,
                                                                       70907, 48812, 71222,
                                                                       28732, 28942, 51062,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 74498, 0, 3,
                                                                       71222, 49037, 71537,
                                                                       28942, 29152, 51377,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 74939, 0, 3,
                                                                       71852, 49487, 72293,
                                                                       29572, 29852, 51692,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 75527, 0, 3,
                                                                       72293, 49802, 72734,
                                                                       29852, 30132, 52112,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 76115, 0, 3,
                                                                       72734, 50117, 73175,
                                                                       30132, 30412, 52532,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 76703, 0, 3,
                                                                       73175, 50432, 73616,
                                                                       30412, 30692, 52952,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 77291, 0, 3,
                                                                       73616, 50747, 74057,
                                                                       30692, 30972, 53372,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 77879, 0, 3,
                                                                       74057, 51062, 74498,
                                                                       30972, 31252, 53792,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 78467, 0, 3,
                                                                       74939, 51692, 75527,
                                                                       31812, 32172, 54212,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 79223, 0, 3,
                                                                       75527, 52112, 76115,
                                                                       32172, 32532, 54752,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 79979, 0, 3,
                                                                       76115, 52532, 76703,
                                                                       32532, 32892, 55292,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 80735, 0, 3,
                                                                       76703, 52952, 77291,
                                                                       32892, 33252, 55832,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 81491, 0, 3,
                                                                       77291, 53372, 77879,
                                                                       33252, 33612, 56372,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 82247, 0, 3,
                                                                       78467, 54212, 79223,
                                                                       34332, 34782, 56912,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 83192, 0, 3,
                                                                       79223, 54752, 79979,
                                                                       34782, 35232, 57587,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 84137, 0, 3,
                                                                       79979, 55292, 80735,
                                                                       35232, 35682, 58262,
                                                                       ncols, gamma, p, q);

                    compute_prim_slh_three_center_electron_repulsion_0(buffer, 85082, 0, 3,
                                                                       80735, 55832, 81491,
                                                                       35682, 36132, 58937,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 86027, 0, 3,
                                                                       82247, 56912, 83192,
                                                                       37032, 37582, 59612,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 87182, 0, 3,
                                                                       83192, 57587, 84137,
                                                                       37582, 38132, 60437,
                                                                       ncols, gamma, p, q);

                    compute_prim_smh_three_center_electron_repulsion_0(buffer, 88337, 0, 3,
                                                                       84137, 58262, 85082,
                                                                       38132, 38682, 61262,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 89492, 0, 3,
                                                                       86027, 59612, 87182,
                                                                       39782, 40442, 62087,
                                                                       ncols, gamma, p, q);

                    compute_prim_snh_three_center_electron_repulsion_0(buffer, 90878, 0, 3,
                                                                       87182, 60437, 88337,
                                                                       40442, 41102, 63077,
                                                                       ncols, gamma, p, q);

                    compute_prim_soh_three_center_electron_repulsion_0(buffer, 92264, 0, 3,
                                                                       89492, 62087, 90878,
                                                                       42422, 43202, 64067,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 93902, 74939, 588, ncols);

                    simdfunc::contract_primitives(buffer, 94798, 78467, 756, ncols);

                    simdfunc::contract_primitives(buffer, 95950, 82247, 945, ncols);

                    simdfunc::contract_primitives(buffer, 97390, 86027, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 99150, 89492, 1386, ncols);

                    simdfunc::contract_primitives(buffer, 101262, 92264, 1638, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 94490, 93902, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 95554, 94798, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 96895, 95950, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 98545, 97390, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 100536, 99150, 66, 1, nmax);

        simdtrf::transform_h_inner(buffer, 102900, 101262, 78, 1, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 103758, 94490, 95554, 11, nmax);

        simdtrf::compute_hrr_pk(buffer, coordinates, 104682, 95554, 96895, 11, nmax);

        simdtrf::compute_hrr_pl(buffer, coordinates, 105870, 96895, 98545, 11, nmax);

        simdtrf::compute_hrr_pm(buffer, coordinates, 107355, 98545, 100536, 11, nmax);

        simdtrf::compute_hrr_pn(buffer, coordinates, 109170, 100536, 102900, 11, nmax);

        simdtrf::compute_hrr_di(buffer, coordinates, 111348, 103758, 104682, 11, nmax);

        simdtrf::compute_hrr_dk(buffer, coordinates, 113196, 104682, 105870, 11, nmax);

        simdtrf::compute_hrr_dl(buffer, coordinates, 115572, 105870, 107355, 11, nmax);

        simdtrf::compute_hrr_dm(buffer, coordinates, 118542, 107355, 109170, 11, nmax);

        simdtrf::compute_hrr_fi(buffer, coordinates, 122172, 111348, 113196, 11, nmax);

        simdtrf::compute_hrr_fk(buffer, coordinates, 125252, 113196, 115572, 11, nmax);

        simdtrf::compute_hrr_fl(buffer, coordinates, 129212, 115572, 118542, 11, nmax);

        simdtrf::compute_hrr_gi(buffer, coordinates, 134162, 122172, 125252, 11, nmax);

        simdtrf::compute_hrr_gk(buffer, coordinates, 138782, 125252, 129212, 11, nmax);

        simdtrf::compute_hrr_hi(buffer, coordinates, 144722, 134162, 138782, 11, nmax);

        simdtrf::transform_i_inner(buffer, 151190, 144722, 21, 11, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 151190, 143, nmax);
    }

    for (size_t m = 0; m < 1573; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
