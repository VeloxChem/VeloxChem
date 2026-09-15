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


#include "SimdThreeCenterElectronRepulsionRecFGL.hpp"

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
#include "SimdTransferDG.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferFG.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_fgl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_fgl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 102874, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1071 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 102874, 85658, 5588, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1298, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1301, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1304, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1307, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1310, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1313, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1316, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1319, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1322, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1325, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1328, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1331, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1334, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1337, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1340, 3, 10, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1349, 3, 11, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1358, 3, 12, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1367, 3, 13, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1376, 3, 14, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1385, 3, 15, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1394, 3, 16, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1403, 3, 17, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1412, 3, 18, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1421, 3, 19, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1430, 3, 20, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1439, 3, 21, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1448, 3, 22, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1457, 3, 30, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1475, 3, 33, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1493, 3, 36, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1511, 3, 39, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1529, 3, 42, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1547, 3, 45, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1565, 3, 48, 117,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1583, 3, 51, 123,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1601, 3, 54, 129,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1619, 3, 57, 135,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1637, 3, 60, 141,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1655, 3, 63, 147,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1673, 3, 81, 173,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1703, 3, 87, 183,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1733, 3, 93, 193,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1763, 3, 99, 203,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1793, 3, 105, 213,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1823, 3, 111, 223,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1853, 3, 117, 233,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1883, 3, 123, 243,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1913, 3, 129, 253,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1943, 3, 135, 263,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1973, 3, 141, 273,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2003, 3, 173, 313,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2048, 3, 183, 328,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2093, 3, 193, 343,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2138, 3, 203, 358,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2183, 3, 213, 373,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2228, 3, 223, 388,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2273, 3, 233, 403,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2318, 3, 243, 418,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2363, 3, 253, 433,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2408, 3, 263, 448,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2453, 3, 313, 505,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2516, 3, 328, 526,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2579, 3, 343, 547,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2642, 3, 358, 568,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2705, 3, 373, 589,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2768, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2831, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2894, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 2957, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3020, 3, 505, 750,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3104, 3, 526, 778,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3188, 3, 547, 806,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3272, 3, 568, 834,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3356, 3, 589, 862,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3440, 3, 610, 890,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3524, 3, 631, 918,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 3608, 3, 652, 946,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3692, 3, 750,
                                                                       1046, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3800, 3, 778,
                                                                       1082, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 3908, 3, 806,
                                                                       1118, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4016, 3, 834,
                                                                       1154, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4124, 3, 862,
                                                                       1190, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4232, 3, 890,
                                                                       1226, ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 4340, 3, 918,
                                                                       1262, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4448, 3, 8, 9,
                                                                       1298, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4454, 3, 9, 10,
                                                                       1301, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4460, 3, 10, 11,
                                                                       1304, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4466, 3, 11, 12,
                                                                       1307, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4472, 3, 12, 13,
                                                                       1310, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4478, 3, 13, 14,
                                                                       1313, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4484, 3, 14, 15,
                                                                       1316, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4490, 3, 15, 16,
                                                                       1319, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4496, 3, 16, 17,
                                                                       1322, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4502, 3, 17, 18,
                                                                       1325, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4508, 3, 18, 19,
                                                                       1328, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4514, 3, 19, 20,
                                                                       1331, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4520, 3, 20, 21,
                                                                       1334, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4526, 3, 21, 22,
                                                                       1337, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4532, 0, 3, 4448,
                                                                       1298, 4454, 1340, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4550, 0, 3, 4454,
                                                                       1301, 4460, 1349, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4568, 0, 3, 4460,
                                                                       1304, 4466, 1358, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4586, 0, 3, 4466,
                                                                       1307, 4472, 1367, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4604, 0, 3, 4472,
                                                                       1310, 4478, 1376, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4622, 0, 3, 4478,
                                                                       1313, 4484, 1385, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4640, 0, 3, 4484,
                                                                       1316, 4490, 1394, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4658, 0, 3, 4490,
                                                                       1319, 4496, 1403, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4676, 0, 3, 4496,
                                                                       1322, 4502, 1412, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4694, 0, 3, 4502,
                                                                       1325, 4508, 1421, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4712, 0, 3, 4508,
                                                                       1328, 4514, 1430, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4730, 0, 3, 4514,
                                                                       1331, 4520, 1439, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4748, 0, 3, 4520,
                                                                       1334, 4526, 1448, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4766, 0, 3, 4532,
                                                                       1340, 4550, 69, 75, 1457,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4802, 0, 3, 4550,
                                                                       1349, 4568, 75, 81, 1475,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4838, 0, 3, 4568,
                                                                       1358, 4586, 81, 87, 1493,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4874, 0, 3, 4586,
                                                                       1367, 4604, 87, 93, 1511,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4910, 0, 3, 4604,
                                                                       1376, 4622, 93, 99, 1529,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4946, 0, 3, 4622,
                                                                       1385, 4640, 99, 105, 1547,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4982, 0, 3, 4640,
                                                                       1394, 4658, 105, 111,
                                                                       1565, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5018, 0, 3, 4658,
                                                                       1403, 4676, 111, 117,
                                                                       1583, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5054, 0, 3, 4676,
                                                                       1412, 4694, 117, 123,
                                                                       1601, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5090, 0, 3, 4694,
                                                                       1421, 4712, 123, 129,
                                                                       1619, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5126, 0, 3, 4712,
                                                                       1430, 4730, 129, 135,
                                                                       1637, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5162, 0, 3, 4730,
                                                                       1439, 4748, 135, 141,
                                                                       1655, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5198, 0, 3, 4766,
                                                                       1457, 4802, 153, 163,
                                                                       1673, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5258, 0, 3, 4802,
                                                                       1475, 4838, 163, 173,
                                                                       1703, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5318, 0, 3, 4838,
                                                                       1493, 4874, 173, 183,
                                                                       1733, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5378, 0, 3, 4874,
                                                                       1511, 4910, 183, 193,
                                                                       1763, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5438, 0, 3, 4910,
                                                                       1529, 4946, 193, 203,
                                                                       1793, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5498, 0, 3, 4946,
                                                                       1547, 4982, 203, 213,
                                                                       1823, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5558, 0, 3, 4982,
                                                                       1565, 5018, 213, 223,
                                                                       1853, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5618, 0, 3, 5018,
                                                                       1583, 5054, 223, 233,
                                                                       1883, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5678, 0, 3, 5054,
                                                                       1601, 5090, 233, 243,
                                                                       1913, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5738, 0, 3, 5090,
                                                                       1619, 5126, 243, 253,
                                                                       1943, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5798, 0, 3, 5126,
                                                                       1637, 5162, 253, 263,
                                                                       1973, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5858, 0, 3, 5198,
                                                                       1673, 5258, 283, 298,
                                                                       2003, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5948, 0, 3, 5258,
                                                                       1703, 5318, 298, 313,
                                                                       2048, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6038, 0, 3, 5318,
                                                                       1733, 5378, 313, 328,
                                                                       2093, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6128, 0, 3, 5378,
                                                                       1763, 5438, 328, 343,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6218, 0, 3, 5438,
                                                                       1793, 5498, 343, 358,
                                                                       2183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6308, 0, 3, 5498,
                                                                       1823, 5558, 358, 373,
                                                                       2228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6398, 0, 3, 5558,
                                                                       1853, 5618, 373, 388,
                                                                       2273, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6488, 0, 3, 5618,
                                                                       1883, 5678, 388, 403,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6578, 0, 3, 5678,
                                                                       1913, 5738, 403, 418,
                                                                       2363, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6668, 0, 3, 5738,
                                                                       1943, 5798, 418, 433,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6758, 0, 3, 5858,
                                                                       2003, 5948, 463, 484,
                                                                       2453, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 6884, 0, 3, 5948,
                                                                       2048, 6038, 484, 505,
                                                                       2516, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7010, 0, 3, 6038,
                                                                       2093, 6128, 505, 526,
                                                                       2579, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7136, 0, 3, 6128,
                                                                       2138, 6218, 526, 547,
                                                                       2642, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7262, 0, 3, 6218,
                                                                       2183, 6308, 547, 568,
                                                                       2705, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7388, 0, 3, 6308,
                                                                       2228, 6398, 568, 589,
                                                                       2768, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7514, 0, 3, 6398,
                                                                       2273, 6488, 589, 610,
                                                                       2831, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7640, 0, 3, 6488,
                                                                       2318, 6578, 610, 631,
                                                                       2894, ncols, gamma, p,
                                                                       q);

                    compute_prim_shd_three_center_electron_repulsion_0(buffer, 7766, 0, 3, 6578,
                                                                       2363, 6668, 631, 652,
                                                                       2957, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 7892, 0, 3, 6758,
                                                                       2453, 6884, 694, 722,
                                                                       3020, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8060, 0, 3, 6884,
                                                                       2516, 7010, 722, 750,
                                                                       3104, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8228, 0, 3, 7010,
                                                                       2579, 7136, 750, 778,
                                                                       3188, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8396, 0, 3, 7136,
                                                                       2642, 7262, 778, 806,
                                                                       3272, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8564, 0, 3, 7262,
                                                                       2705, 7388, 806, 834,
                                                                       3356, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8732, 0, 3, 7388,
                                                                       2768, 7514, 834, 862,
                                                                       3440, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 8900, 0, 3, 7514,
                                                                       2831, 7640, 862, 890,
                                                                       3524, ncols, gamma, p,
                                                                       q);

                    compute_prim_sid_three_center_electron_repulsion_0(buffer, 9068, 0, 3, 7640,
                                                                       2894, 7766, 890, 918,
                                                                       3608, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9236, 0, 3, 7892,
                                                                       3020, 8060, 974, 1010,
                                                                       3692, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9452, 0, 3, 8060,
                                                                       3104, 8228, 1010, 1046,
                                                                       3800, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9668, 0, 3, 8228,
                                                                       3188, 8396, 1046, 1082,
                                                                       3908, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 9884, 0, 3, 8396,
                                                                       3272, 8564, 1082, 1118,
                                                                       4016, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10100, 0, 3, 8564,
                                                                       3356, 8732, 1118, 1154,
                                                                       4124, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10316, 0, 3, 8732,
                                                                       3440, 8900, 1154, 1190,
                                                                       4232, ncols, gamma, p,
                                                                       q);

                    compute_prim_skd_three_center_electron_repulsion_0(buffer, 10532, 0, 3, 8900,
                                                                       3524, 9068, 1190, 1226,
                                                                       4340, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10748, 3, 1298,
                                                                       1301, 4460, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10758, 3, 1301,
                                                                       1304, 4466, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10768, 3, 1304,
                                                                       1307, 4472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10778, 3, 1307,
                                                                       1310, 4478, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10788, 3, 1310,
                                                                       1313, 4484, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10798, 3, 1313,
                                                                       1316, 4490, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10808, 3, 1316,
                                                                       1319, 4496, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10818, 3, 1319,
                                                                       1322, 4502, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10828, 3, 1322,
                                                                       1325, 4508, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10838, 3, 1325,
                                                                       1328, 4514, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10848, 3, 1328,
                                                                       1331, 4520, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10858, 3, 1331,
                                                                       1334, 4526, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10868, 0, 3,
                                                                       10748, 4460, 10758, 4568,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10898, 0, 3,
                                                                       10758, 4466, 10768, 4586,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10928, 0, 3,
                                                                       10768, 4472, 10778, 4604,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10958, 0, 3,
                                                                       10778, 4478, 10788, 4622,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10988, 0, 3,
                                                                       10788, 4484, 10798, 4640,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11018, 0, 3,
                                                                       10798, 4490, 10808, 4658,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11048, 0, 3,
                                                                       10808, 4496, 10818, 4676,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11078, 0, 3,
                                                                       10818, 4502, 10828, 4694,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11108, 0, 3,
                                                                       10828, 4508, 10838, 4712,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11138, 0, 3,
                                                                       10838, 4514, 10848, 4730,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 11168, 0, 3,
                                                                       10848, 4520, 10858, 4748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11198, 0, 3,
                                                                       10868, 4568, 10898, 1457,
                                                                       1475, 4838, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11258, 0, 3,
                                                                       10898, 4586, 10928, 1475,
                                                                       1493, 4874, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11318, 0, 3,
                                                                       10928, 4604, 10958, 1493,
                                                                       1511, 4910, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11378, 0, 3,
                                                                       10958, 4622, 10988, 1511,
                                                                       1529, 4946, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11438, 0, 3,
                                                                       10988, 4640, 11018, 1529,
                                                                       1547, 4982, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11498, 0, 3,
                                                                       11018, 4658, 11048, 1547,
                                                                       1565, 5018, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11558, 0, 3,
                                                                       11048, 4676, 11078, 1565,
                                                                       1583, 5054, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11618, 0, 3,
                                                                       11078, 4694, 11108, 1583,
                                                                       1601, 5090, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11678, 0, 3,
                                                                       11108, 4712, 11138, 1601,
                                                                       1619, 5126, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 11738, 0, 3,
                                                                       11138, 4730, 11168, 1619,
                                                                       1637, 5162, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11798, 0, 3,
                                                                       11198, 4838, 11258, 1673,
                                                                       1703, 5318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11898, 0, 3,
                                                                       11258, 4874, 11318, 1703,
                                                                       1733, 5378, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11998, 0, 3,
                                                                       11318, 4910, 11378, 1733,
                                                                       1763, 5438, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12098, 0, 3,
                                                                       11378, 4946, 11438, 1763,
                                                                       1793, 5498, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12198, 0, 3,
                                                                       11438, 4982, 11498, 1793,
                                                                       1823, 5558, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12298, 0, 3,
                                                                       11498, 5018, 11558, 1823,
                                                                       1853, 5618, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12398, 0, 3,
                                                                       11558, 5054, 11618, 1853,
                                                                       1883, 5678, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12498, 0, 3,
                                                                       11618, 5090, 11678, 1883,
                                                                       1913, 5738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 12598, 0, 3,
                                                                       11678, 5126, 11738, 1913,
                                                                       1943, 5798, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12698, 0, 3,
                                                                       11798, 5318, 11898, 2003,
                                                                       2048, 6038, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12848, 0, 3,
                                                                       11898, 5378, 11998, 2048,
                                                                       2093, 6128, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 12998, 0, 3,
                                                                       11998, 5438, 12098, 2093,
                                                                       2138, 6218, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13148, 0, 3,
                                                                       12098, 5498, 12198, 2138,
                                                                       2183, 6308, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13298, 0, 3,
                                                                       12198, 5558, 12298, 2183,
                                                                       2228, 6398, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13448, 0, 3,
                                                                       12298, 5618, 12398, 2228,
                                                                       2273, 6488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13598, 0, 3,
                                                                       12398, 5678, 12498, 2273,
                                                                       2318, 6578, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13748, 0, 3,
                                                                       12498, 5738, 12598, 2318,
                                                                       2363, 6668, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 13898, 0, 3,
                                                                       12698, 6038, 12848, 2453,
                                                                       2516, 7010, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14108, 0, 3,
                                                                       12848, 6128, 12998, 2516,
                                                                       2579, 7136, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14318, 0, 3,
                                                                       12998, 6218, 13148, 2579,
                                                                       2642, 7262, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14528, 0, 3,
                                                                       13148, 6308, 13298, 2642,
                                                                       2705, 7388, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14738, 0, 3,
                                                                       13298, 6398, 13448, 2705,
                                                                       2768, 7514, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 14948, 0, 3,
                                                                       13448, 6488, 13598, 2768,
                                                                       2831, 7640, ncols, gamma,
                                                                       p, q);

                    compute_prim_shf_three_center_electron_repulsion_0(buffer, 15158, 0, 3,
                                                                       13598, 6578, 13748, 2831,
                                                                       2894, 7766, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 15368, 0, 3,
                                                                       13898, 7010, 14108, 3020,
                                                                       3104, 8228, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 15648, 0, 3,
                                                                       14108, 7136, 14318, 3104,
                                                                       3188, 8396, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 15928, 0, 3,
                                                                       14318, 7262, 14528, 3188,
                                                                       3272, 8564, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16208, 0, 3,
                                                                       14528, 7388, 14738, 3272,
                                                                       3356, 8732, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16488, 0, 3,
                                                                       14738, 7514, 14948, 3356,
                                                                       3440, 8900, ncols, gamma,
                                                                       p, q);

                    compute_prim_sif_three_center_electron_repulsion_0(buffer, 16768, 0, 3,
                                                                       14948, 7640, 15158, 3440,
                                                                       3524, 9068, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 17048, 0, 3,
                                                                       15368, 8228, 15648, 3692,
                                                                       3800, 9668, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 17408, 0, 3,
                                                                       15648, 8396, 15928, 3800,
                                                                       3908, 9884, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 17768, 0, 3,
                                                                       15928, 8564, 16208, 3908,
                                                                       4016, 10100, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 18128, 0, 3,
                                                                       16208, 8732, 16488, 4016,
                                                                       4124, 10316, ncols, gamma,
                                                                       p, q);

                    compute_prim_skf_three_center_electron_repulsion_0(buffer, 18488, 0, 3,
                                                                       16488, 8900, 16768, 4124,
                                                                       4232, 10532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18848, 3, 4448,
                                                                       4454, 10748, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18863, 3, 4454,
                                                                       4460, 10758, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18878, 3, 4460,
                                                                       4466, 10768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18893, 3, 4466,
                                                                       4472, 10778, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18908, 3, 4472,
                                                                       4478, 10788, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18923, 3, 4478,
                                                                       4484, 10798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18938, 3, 4484,
                                                                       4490, 10808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18953, 3, 4490,
                                                                       4496, 10818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18968, 3, 4496,
                                                                       4502, 10828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18983, 3, 4502,
                                                                       4508, 10838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18998, 3, 4508,
                                                                       4514, 10848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 19013, 3, 4514,
                                                                       4520, 10858, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19028, 0, 3,
                                                                       18848, 10748, 18863, 4532,
                                                                       4550, 10868, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19073, 0, 3,
                                                                       18863, 10758, 18878, 4550,
                                                                       4568, 10898, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19118, 0, 3,
                                                                       18878, 10768, 18893, 4568,
                                                                       4586, 10928, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19163, 0, 3,
                                                                       18893, 10778, 18908, 4586,
                                                                       4604, 10958, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19208, 0, 3,
                                                                       18908, 10788, 18923, 4604,
                                                                       4622, 10988, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19253, 0, 3,
                                                                       18923, 10798, 18938, 4622,
                                                                       4640, 11018, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19298, 0, 3,
                                                                       18938, 10808, 18953, 4640,
                                                                       4658, 11048, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19343, 0, 3,
                                                                       18953, 10818, 18968, 4658,
                                                                       4676, 11078, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19388, 0, 3,
                                                                       18968, 10828, 18983, 4676,
                                                                       4694, 11108, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19433, 0, 3,
                                                                       18983, 10838, 18998, 4694,
                                                                       4712, 11138, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19478, 0, 3,
                                                                       18998, 10848, 19013, 4712,
                                                                       4730, 11168, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19523, 0, 3,
                                                                       19028, 10868, 19073, 4766,
                                                                       4802, 11198, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19613, 0, 3,
                                                                       19073, 10898, 19118, 4802,
                                                                       4838, 11258, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19703, 0, 3,
                                                                       19118, 10928, 19163, 4838,
                                                                       4874, 11318, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19793, 0, 3,
                                                                       19163, 10958, 19208, 4874,
                                                                       4910, 11378, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19883, 0, 3,
                                                                       19208, 10988, 19253, 4910,
                                                                       4946, 11438, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19973, 0, 3,
                                                                       19253, 11018, 19298, 4946,
                                                                       4982, 11498, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20063, 0, 3,
                                                                       19298, 11048, 19343, 4982,
                                                                       5018, 11558, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20153, 0, 3,
                                                                       19343, 11078, 19388, 5018,
                                                                       5054, 11618, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20243, 0, 3,
                                                                       19388, 11108, 19433, 5054,
                                                                       5090, 11678, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20333, 0, 3,
                                                                       19433, 11138, 19478, 5090,
                                                                       5126, 11738, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20423, 0, 3,
                                                                       19523, 11198, 19613, 5198,
                                                                       5258, 11798, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20573, 0, 3,
                                                                       19613, 11258, 19703, 5258,
                                                                       5318, 11898, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20723, 0, 3,
                                                                       19703, 11318, 19793, 5318,
                                                                       5378, 11998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20873, 0, 3,
                                                                       19793, 11378, 19883, 5378,
                                                                       5438, 12098, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21023, 0, 3,
                                                                       19883, 11438, 19973, 5438,
                                                                       5498, 12198, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21173, 0, 3,
                                                                       19973, 11498, 20063, 5498,
                                                                       5558, 12298, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21323, 0, 3,
                                                                       20063, 11558, 20153, 5558,
                                                                       5618, 12398, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21473, 0, 3,
                                                                       20153, 11618, 20243, 5618,
                                                                       5678, 12498, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 21623, 0, 3,
                                                                       20243, 11678, 20333, 5678,
                                                                       5738, 12598, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21773, 0, 3,
                                                                       20423, 11798, 20573, 5858,
                                                                       5948, 12698, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21998, 0, 3,
                                                                       20573, 11898, 20723, 5948,
                                                                       6038, 12848, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22223, 0, 3,
                                                                       20723, 11998, 20873, 6038,
                                                                       6128, 12998, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22448, 0, 3,
                                                                       20873, 12098, 21023, 6128,
                                                                       6218, 13148, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22673, 0, 3,
                                                                       21023, 12198, 21173, 6218,
                                                                       6308, 13298, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 22898, 0, 3,
                                                                       21173, 12298, 21323, 6308,
                                                                       6398, 13448, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23123, 0, 3,
                                                                       21323, 12398, 21473, 6398,
                                                                       6488, 13598, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 23348, 0, 3,
                                                                       21473, 12498, 21623, 6488,
                                                                       6578, 13748, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 23573, 0, 3,
                                                                       21773, 12698, 21998, 6758,
                                                                       6884, 13898, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 23888, 0, 3,
                                                                       21998, 12848, 22223, 6884,
                                                                       7010, 14108, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24203, 0, 3,
                                                                       22223, 12998, 22448, 7010,
                                                                       7136, 14318, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24518, 0, 3,
                                                                       22448, 13148, 22673, 7136,
                                                                       7262, 14528, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 24833, 0, 3,
                                                                       22673, 13298, 22898, 7262,
                                                                       7388, 14738, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25148, 0, 3,
                                                                       22898, 13448, 23123, 7388,
                                                                       7514, 14948, ncols, gamma,
                                                                       p, q);

                    compute_prim_shg_three_center_electron_repulsion_0(buffer, 25463, 0, 3,
                                                                       23123, 13598, 23348, 7514,
                                                                       7640, 15158, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 25778, 0, 3,
                                                                       23573, 13898, 23888, 7892,
                                                                       8060, 15368, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 26198, 0, 3,
                                                                       23888, 14108, 24203, 8060,
                                                                       8228, 15648, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 26618, 0, 3,
                                                                       24203, 14318, 24518, 8228,
                                                                       8396, 15928, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 27038, 0, 3,
                                                                       24518, 14528, 24833, 8396,
                                                                       8564, 16208, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 27458, 0, 3,
                                                                       24833, 14738, 25148, 8564,
                                                                       8732, 16488, ncols, gamma,
                                                                       p, q);

                    compute_prim_sig_three_center_electron_repulsion_0(buffer, 27878, 0, 3,
                                                                       25148, 14948, 25463, 8732,
                                                                       8900, 16768, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 28298, 0, 3,
                                                                       25778, 15368, 26198, 9236,
                                                                       9452, 17048, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 28838, 0, 3,
                                                                       26198, 15648, 26618, 9452,
                                                                       9668, 17408, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 29378, 0, 3,
                                                                       26618, 15928, 27038, 9668,
                                                                       9884, 17768, ncols, gamma,
                                                                       p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 29918, 0, 3,
                                                                       27038, 16208, 27458, 9884,
                                                                       10100, 18128, ncols,
                                                                       gamma, p, q);

                    compute_prim_skg_three_center_electron_repulsion_0(buffer, 30458, 0, 3,
                                                                       27458, 16488, 27878,
                                                                       10100, 10316, 18488,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 30998, 3, 10748,
                                                                       10758, 18878, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31019, 3, 10758,
                                                                       10768, 18893, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31040, 3, 10768,
                                                                       10778, 18908, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31061, 3, 10778,
                                                                       10788, 18923, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31082, 3, 10788,
                                                                       10798, 18938, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31103, 3, 10798,
                                                                       10808, 18953, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31124, 3, 10808,
                                                                       10818, 18968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31145, 3, 10818,
                                                                       10828, 18983, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31166, 3, 10828,
                                                                       10838, 18998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 31187, 3, 10838,
                                                                       10848, 19013, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31208, 0, 3,
                                                                       30998, 18878, 31019,
                                                                       10868, 10898, 19118,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31271, 0, 3,
                                                                       31019, 18893, 31040,
                                                                       10898, 10928, 19163,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31334, 0, 3,
                                                                       31040, 18908, 31061,
                                                                       10928, 10958, 19208,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31397, 0, 3,
                                                                       31061, 18923, 31082,
                                                                       10958, 10988, 19253,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31460, 0, 3,
                                                                       31082, 18938, 31103,
                                                                       10988, 11018, 19298,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31523, 0, 3,
                                                                       31103, 18953, 31124,
                                                                       11018, 11048, 19343,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31586, 0, 3,
                                                                       31124, 18968, 31145,
                                                                       11048, 11078, 19388,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31649, 0, 3,
                                                                       31145, 18983, 31166,
                                                                       11078, 11108, 19433,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 31712, 0, 3,
                                                                       31166, 18998, 31187,
                                                                       11108, 11138, 19478,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31775, 0, 3,
                                                                       31208, 19118, 31271,
                                                                       11198, 11258, 19703,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 31901, 0, 3,
                                                                       31271, 19163, 31334,
                                                                       11258, 11318, 19793,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32027, 0, 3,
                                                                       31334, 19208, 31397,
                                                                       11318, 11378, 19883,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32153, 0, 3,
                                                                       31397, 19253, 31460,
                                                                       11378, 11438, 19973,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32279, 0, 3,
                                                                       31460, 19298, 31523,
                                                                       11438, 11498, 20063,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32405, 0, 3,
                                                                       31523, 19343, 31586,
                                                                       11498, 11558, 20153,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32531, 0, 3,
                                                                       31586, 19388, 31649,
                                                                       11558, 11618, 20243,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 32657, 0, 3,
                                                                       31649, 19433, 31712,
                                                                       11618, 11678, 20333,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 32783, 0, 3,
                                                                       31775, 19703, 31901,
                                                                       11798, 11898, 20723,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 32993, 0, 3,
                                                                       31901, 19793, 32027,
                                                                       11898, 11998, 20873,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33203, 0, 3,
                                                                       32027, 19883, 32153,
                                                                       11998, 12098, 21023,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33413, 0, 3,
                                                                       32153, 19973, 32279,
                                                                       12098, 12198, 21173,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33623, 0, 3,
                                                                       32279, 20063, 32405,
                                                                       12198, 12298, 21323,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 33833, 0, 3,
                                                                       32405, 20153, 32531,
                                                                       12298, 12398, 21473,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 34043, 0, 3,
                                                                       32531, 20243, 32657,
                                                                       12398, 12498, 21623,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 34253, 0, 3,
                                                                       32783, 20723, 32993,
                                                                       12698, 12848, 22223,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 34568, 0, 3,
                                                                       32993, 20873, 33203,
                                                                       12848, 12998, 22448,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 34883, 0, 3,
                                                                       33203, 21023, 33413,
                                                                       12998, 13148, 22673,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 35198, 0, 3,
                                                                       33413, 21173, 33623,
                                                                       13148, 13298, 22898,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 35513, 0, 3,
                                                                       33623, 21323, 33833,
                                                                       13298, 13448, 23123,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 35828, 0, 3,
                                                                       33833, 21473, 34043,
                                                                       13448, 13598, 23348,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 36143, 0, 3,
                                                                       34253, 22223, 34568,
                                                                       13898, 14108, 24203,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 36584, 0, 3,
                                                                       34568, 22448, 34883,
                                                                       14108, 14318, 24518,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 37025, 0, 3,
                                                                       34883, 22673, 35198,
                                                                       14318, 14528, 24833,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 37466, 0, 3,
                                                                       35198, 22898, 35513,
                                                                       14528, 14738, 25148,
                                                                       ncols, gamma, p, q);

                    compute_prim_shh_three_center_electron_repulsion_0(buffer, 37907, 0, 3,
                                                                       35513, 23123, 35828,
                                                                       14738, 14948, 25463,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 38348, 0, 3,
                                                                       36143, 24203, 36584,
                                                                       15368, 15648, 26618,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 38936, 0, 3,
                                                                       36584, 24518, 37025,
                                                                       15648, 15928, 27038,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 39524, 0, 3,
                                                                       37025, 24833, 37466,
                                                                       15928, 16208, 27458,
                                                                       ncols, gamma, p, q);

                    compute_prim_sih_three_center_electron_repulsion_0(buffer, 40112, 0, 3,
                                                                       37466, 25148, 37907,
                                                                       16208, 16488, 27878,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 40700, 0, 3,
                                                                       38348, 26618, 38936,
                                                                       17048, 17408, 29378,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 41456, 0, 3,
                                                                       38936, 27038, 39524,
                                                                       17408, 17768, 29918,
                                                                       ncols, gamma, p, q);

                    compute_prim_skh_three_center_electron_repulsion_0(buffer, 42212, 0, 3,
                                                                       39524, 27458, 40112,
                                                                       17768, 18128, 30458,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42968, 3, 18848,
                                                                       18863, 30998, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42996, 3, 18863,
                                                                       18878, 31019, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43024, 3, 18878,
                                                                       18893, 31040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43052, 3, 18893,
                                                                       18908, 31061, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43080, 3, 18908,
                                                                       18923, 31082, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43108, 3, 18923,
                                                                       18938, 31103, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43136, 3, 18938,
                                                                       18953, 31124, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43164, 3, 18953,
                                                                       18968, 31145, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43192, 3, 18968,
                                                                       18983, 31166, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43220, 3, 18983,
                                                                       18998, 31187, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43248, 0, 3,
                                                                       42968, 30998, 42996,
                                                                       19028, 19073, 31208,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43332, 0, 3,
                                                                       42996, 31019, 43024,
                                                                       19073, 19118, 31271,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43416, 0, 3,
                                                                       43024, 31040, 43052,
                                                                       19118, 19163, 31334,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43500, 0, 3,
                                                                       43052, 31061, 43080,
                                                                       19163, 19208, 31397,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43584, 0, 3,
                                                                       43080, 31082, 43108,
                                                                       19208, 19253, 31460,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43668, 0, 3,
                                                                       43108, 31103, 43136,
                                                                       19253, 19298, 31523,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43752, 0, 3,
                                                                       43136, 31124, 43164,
                                                                       19298, 19343, 31586,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43836, 0, 3,
                                                                       43164, 31145, 43192,
                                                                       19343, 19388, 31649,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 43920, 0, 3,
                                                                       43192, 31166, 43220,
                                                                       19388, 19433, 31712,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44004, 0, 3,
                                                                       43248, 31208, 43332,
                                                                       19523, 19613, 31775,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44172, 0, 3,
                                                                       43332, 31271, 43416,
                                                                       19613, 19703, 31901,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44340, 0, 3,
                                                                       43416, 31334, 43500,
                                                                       19703, 19793, 32027,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44508, 0, 3,
                                                                       43500, 31397, 43584,
                                                                       19793, 19883, 32153,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44676, 0, 3,
                                                                       43584, 31460, 43668,
                                                                       19883, 19973, 32279,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 44844, 0, 3,
                                                                       43668, 31523, 43752,
                                                                       19973, 20063, 32405,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 45012, 0, 3,
                                                                       43752, 31586, 43836,
                                                                       20063, 20153, 32531,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 45180, 0, 3,
                                                                       43836, 31649, 43920,
                                                                       20153, 20243, 32657,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 45348, 0, 3,
                                                                       44004, 31775, 44172,
                                                                       20423, 20573, 32783,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 45628, 0, 3,
                                                                       44172, 31901, 44340,
                                                                       20573, 20723, 32993,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 45908, 0, 3,
                                                                       44340, 32027, 44508,
                                                                       20723, 20873, 33203,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 46188, 0, 3,
                                                                       44508, 32153, 44676,
                                                                       20873, 21023, 33413,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 46468, 0, 3,
                                                                       44676, 32279, 44844,
                                                                       21023, 21173, 33623,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 46748, 0, 3,
                                                                       44844, 32405, 45012,
                                                                       21173, 21323, 33833,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 47028, 0, 3,
                                                                       45012, 32531, 45180,
                                                                       21323, 21473, 34043,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 47308, 0, 3,
                                                                       45348, 32783, 45628,
                                                                       21773, 21998, 34253,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 47728, 0, 3,
                                                                       45628, 32993, 45908,
                                                                       21998, 22223, 34568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 48148, 0, 3,
                                                                       45908, 33203, 46188,
                                                                       22223, 22448, 34883,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 48568, 0, 3,
                                                                       46188, 33413, 46468,
                                                                       22448, 22673, 35198,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 48988, 0, 3,
                                                                       46468, 33623, 46748,
                                                                       22673, 22898, 35513,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgi_three_center_electron_repulsion_0(buffer, 49408, 0, 3,
                                                                       46748, 33833, 47028,
                                                                       22898, 23123, 35828,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 49828, 0, 3,
                                                                       47308, 34253, 47728,
                                                                       23573, 23888, 36143,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 50416, 0, 3,
                                                                       47728, 34568, 48148,
                                                                       23888, 24203, 36584,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 51004, 0, 3,
                                                                       48148, 34883, 48568,
                                                                       24203, 24518, 37025,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 51592, 0, 3,
                                                                       48568, 35198, 48988,
                                                                       24518, 24833, 37466,
                                                                       ncols, gamma, p, q);

                    compute_prim_shi_three_center_electron_repulsion_0(buffer, 52180, 0, 3,
                                                                       48988, 35513, 49408,
                                                                       24833, 25148, 37907,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 52768, 0, 3,
                                                                       49828, 36143, 50416,
                                                                       25778, 26198, 38348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 53552, 0, 3,
                                                                       50416, 36584, 51004,
                                                                       26198, 26618, 38936,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 54336, 0, 3,
                                                                       51004, 37025, 51592,
                                                                       26618, 27038, 39524,
                                                                       ncols, gamma, p, q);

                    compute_prim_sii_three_center_electron_repulsion_0(buffer, 55120, 0, 3,
                                                                       51592, 37466, 52180,
                                                                       27038, 27458, 40112,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 55904, 0, 3,
                                                                       52768, 38348, 53552,
                                                                       28298, 28838, 40700,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 56912, 0, 3,
                                                                       53552, 38936, 54336,
                                                                       28838, 29378, 41456,
                                                                       ncols, gamma, p, q);

                    compute_prim_ski_three_center_electron_repulsion_0(buffer, 57920, 0, 3,
                                                                       54336, 39524, 55120,
                                                                       29378, 29918, 42212,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 58928, 3, 30998,
                                                                       31019, 43024, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 58964, 3, 31019,
                                                                       31040, 43052, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59000, 3, 31040,
                                                                       31061, 43080, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59036, 3, 31061,
                                                                       31082, 43108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59072, 3, 31082,
                                                                       31103, 43136, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59108, 3, 31103,
                                                                       31124, 43164, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59144, 3, 31124,
                                                                       31145, 43192, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 59180, 3, 31145,
                                                                       31166, 43220, ncols,
                                                                       gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59216, 0, 3,
                                                                       58928, 43024, 58964,
                                                                       31208, 31271, 43416,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59324, 0, 3,
                                                                       58964, 43052, 59000,
                                                                       31271, 31334, 43500,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59432, 0, 3,
                                                                       59000, 43080, 59036,
                                                                       31334, 31397, 43584,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59540, 0, 3,
                                                                       59036, 43108, 59072,
                                                                       31397, 31460, 43668,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59648, 0, 3,
                                                                       59072, 43136, 59108,
                                                                       31460, 31523, 43752,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59756, 0, 3,
                                                                       59108, 43164, 59144,
                                                                       31523, 31586, 43836,
                                                                       ncols, gamma, p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 59864, 0, 3,
                                                                       59144, 43192, 59180,
                                                                       31586, 31649, 43920,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 59972, 0, 3,
                                                                       59216, 43416, 59324,
                                                                       31775, 31901, 44340,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 60188, 0, 3,
                                                                       59324, 43500, 59432,
                                                                       31901, 32027, 44508,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 60404, 0, 3,
                                                                       59432, 43584, 59540,
                                                                       32027, 32153, 44676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 60620, 0, 3,
                                                                       59540, 43668, 59648,
                                                                       32153, 32279, 44844,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 60836, 0, 3,
                                                                       59648, 43752, 59756,
                                                                       32279, 32405, 45012,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 61052, 0, 3,
                                                                       59756, 43836, 59864,
                                                                       32405, 32531, 45180,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 61268, 0, 3,
                                                                       59972, 44340, 60188,
                                                                       32783, 32993, 45908,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 61628, 0, 3,
                                                                       60188, 44508, 60404,
                                                                       32993, 33203, 46188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 61988, 0, 3,
                                                                       60404, 44676, 60620,
                                                                       33203, 33413, 46468,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 62348, 0, 3,
                                                                       60620, 44844, 60836,
                                                                       33413, 33623, 46748,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfk_three_center_electron_repulsion_0(buffer, 62708, 0, 3,
                                                                       60836, 45012, 61052,
                                                                       33623, 33833, 47028,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 63068, 0, 3,
                                                                       61268, 45908, 61628,
                                                                       34253, 34568, 48148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 63608, 0, 3,
                                                                       61628, 46188, 61988,
                                                                       34568, 34883, 48568,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 64148, 0, 3,
                                                                       61988, 46468, 62348,
                                                                       34883, 35198, 48988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgk_three_center_electron_repulsion_0(buffer, 64688, 0, 3,
                                                                       62348, 46748, 62708,
                                                                       35198, 35513, 49408,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 65228, 0, 3,
                                                                       63068, 48148, 63608,
                                                                       36143, 36584, 51004,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 65984, 0, 3,
                                                                       63608, 48568, 64148,
                                                                       36584, 37025, 51592,
                                                                       ncols, gamma, p, q);

                    compute_prim_shk_three_center_electron_repulsion_0(buffer, 66740, 0, 3,
                                                                       64148, 48988, 64688,
                                                                       37025, 37466, 52180,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 67496, 0, 3,
                                                                       65228, 51004, 65984,
                                                                       38348, 38936, 54336,
                                                                       ncols, gamma, p, q);

                    compute_prim_sik_three_center_electron_repulsion_0(buffer, 68504, 0, 3,
                                                                       65984, 51592, 66740,
                                                                       38936, 39524, 55120,
                                                                       ncols, gamma, p, q);

                    compute_prim_skk_three_center_electron_repulsion_0(buffer, 69512, 0, 3,
                                                                       67496, 54336, 68504,
                                                                       40700, 41456, 57920,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 70808, 3, 42968,
                                                                       42996, 58928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 70853, 3, 42996,
                                                                       43024, 58964, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 70898, 3, 43024,
                                                                       43052, 59000, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 70943, 3, 43052,
                                                                       43080, 59036, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 70988, 3, 43080,
                                                                       43108, 59072, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 71033, 3, 43108,
                                                                       43136, 59108, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 71078, 3, 43136,
                                                                       43164, 59144, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 71123, 3, 43164,
                                                                       43192, 59180, ncols,
                                                                       gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71168, 0, 3,
                                                                       70808, 58928, 70853,
                                                                       43248, 43332, 59216,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71303, 0, 3,
                                                                       70853, 58964, 70898,
                                                                       43332, 43416, 59324,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71438, 0, 3,
                                                                       70898, 59000, 70943,
                                                                       43416, 43500, 59432,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71573, 0, 3,
                                                                       70943, 59036, 70988,
                                                                       43500, 43584, 59540,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71708, 0, 3,
                                                                       70988, 59072, 71033,
                                                                       43584, 43668, 59648,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71843, 0, 3,
                                                                       71033, 59108, 71078,
                                                                       43668, 43752, 59756,
                                                                       ncols, gamma, p, q);

                    compute_prim_spl_three_center_electron_repulsion_0(buffer, 71978, 0, 3,
                                                                       71078, 59144, 71123,
                                                                       43752, 43836, 59864,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 72113, 0, 3,
                                                                       71168, 59216, 71303,
                                                                       44004, 44172, 59972,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 72383, 0, 3,
                                                                       71303, 59324, 71438,
                                                                       44172, 44340, 60188,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 72653, 0, 3,
                                                                       71438, 59432, 71573,
                                                                       44340, 44508, 60404,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 72923, 0, 3,
                                                                       71573, 59540, 71708,
                                                                       44508, 44676, 60620,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 73193, 0, 3,
                                                                       71708, 59648, 71843,
                                                                       44676, 44844, 60836,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdl_three_center_electron_repulsion_0(buffer, 73463, 0, 3,
                                                                       71843, 59756, 71978,
                                                                       44844, 45012, 61052,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 73733, 0, 3,
                                                                       72113, 59972, 72383,
                                                                       45348, 45628, 61268,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 74183, 0, 3,
                                                                       72383, 60188, 72653,
                                                                       45628, 45908, 61628,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 74633, 0, 3,
                                                                       72653, 60404, 72923,
                                                                       45908, 46188, 61988,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 75083, 0, 3,
                                                                       72923, 60620, 73193,
                                                                       46188, 46468, 62348,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfl_three_center_electron_repulsion_0(buffer, 75533, 0, 3,
                                                                       73193, 60836, 73463,
                                                                       46468, 46748, 62708,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 75983, 0, 3,
                                                                       73733, 61268, 74183,
                                                                       47308, 47728, 63068,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 76658, 0, 3,
                                                                       74183, 61628, 74633,
                                                                       47728, 48148, 63608,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 77333, 0, 3,
                                                                       74633, 61988, 75083,
                                                                       48148, 48568, 64148,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgl_three_center_electron_repulsion_0(buffer, 78008, 0, 3,
                                                                       75083, 62348, 75533,
                                                                       48568, 48988, 64688,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 78683, 0, 3,
                                                                       75983, 63068, 76658,
                                                                       49828, 50416, 65228,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 79628, 0, 3,
                                                                       76658, 63608, 77333,
                                                                       50416, 51004, 65984,
                                                                       ncols, gamma, p, q);

                    compute_prim_shl_three_center_electron_repulsion_0(buffer, 80573, 0, 3,
                                                                       77333, 64148, 78008,
                                                                       51004, 51592, 66740,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 81518, 0, 3,
                                                                       78683, 65228, 79628,
                                                                       52768, 53552, 67496,
                                                                       ncols, gamma, p, q);

                    compute_prim_sil_three_center_electron_repulsion_0(buffer, 82778, 0, 3,
                                                                       79628, 65984, 80573,
                                                                       53552, 54336, 68504,
                                                                       ncols, gamma, p, q);

                    compute_prim_skl_three_center_electron_repulsion_0(buffer, 84038, 0, 3,
                                                                       81518, 67496, 82778,
                                                                       55904, 56912, 69512,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 85658, 75983, 675, ncols);

                    simdfunc::contract_primitives(buffer, 86588, 78683, 945, ncols);

                    simdfunc::contract_primitives(buffer, 87890, 81518, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 89626, 84038, 1620, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 86333, 85658, 15, 1, nmax);

        simdtrf::transform_l_inner(buffer, 87533, 86588, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 89150, 87890, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 91246, 89626, 36, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 91858, 86333, 87533, 17, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 92623, 87533, 89150, 17, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 93694, 89150, 91246, 17, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 95122, 91858, 92623, 17, nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 96652, 92623, 93694, 17, nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 98794, 95122, 96652, 17, nmax);

        simdtrf::transform_g_inner(buffer, 101344, 98794, 10, 17, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 101344, 153, nmax);
    }

    for (size_t m = 0; m < 1071; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
